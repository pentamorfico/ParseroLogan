# download.py
import json
import subprocess
import threading
from pathlib import Path

from .models import ProcessedSRA
from .annotation import process_annotations
from .hmm_search import search_hmm_streaming
from .output import save_results
from .fasta import count_fasta_contigs

# Create a global variable to hold the manager.
global_manager = None

def set_manager(mgr):
    """Set the global manager. This must be called once in the main process before spawning workers."""
    global global_manager
    global_manager = mgr

def download_sra_to_variable(sra_id: str, updates: dict, done_event, hmm_file: Path, min_length: int, evalue: float, output_dir: Path, fix_circles_flag: bool = False, save_all: bool = False):
    # Use the global manager to create a proxy dict for this SRA.
    updates[sra_id] = global_manager.dict({
        "name": sra_id,
        "status": "Downloading...",
        "downloaded_symbol": "[yellow]?[/yellow]",
        "annotated_symbol": "",
        "searched_symbol": "",
        "data": b"",
        "num_contigs": None,
        "num_proteins": None,
        "num_hits": None
    })
    bucket_key = f"s3://logan-pub/c/{sra_id}/{sra_id}.contigs.fa.zst"
    cmd = f"aws s3 cp {bucket_key} - --no-sign-request --no-progress | seqkit seq -m {min_length}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        data = process.stdout.read()
        updates[sra_id]["data"] = data
    except Exception:
        process.kill()
        updates[sra_id]["status"] = "Error"
        updates[sra_id]["downloaded_symbol"] = "[red]✗[/red]"
        done_event.set()
        return

    ret = process.wait()

    if ret != 0:
        updates[sra_id]["status"] = "Error"
        updates[sra_id]["downloaded_symbol"] = "[red]✗[/red]"
        done_event.set()
        return

    updates[sra_id]["downloaded_symbol"] = "[green]✓[/green]"
    num_contigs = count_fasta_contigs(updates[sra_id]["data"])
    updates[sra_id]["num_contigs"] = num_contigs if num_contigs > 0 else None
    updates[sra_id]["status"] = "Annotating..."

    try:
        all_contigs, all_annotations, all_proteins, hmmer_sequences = process_annotations(
            updates[sra_id]["data"], min_length, sra_id, fix_circles_flag=fix_circles_flag
        )
    except Exception:
        import traceback
        traceback.print_exc()
        updates[sra_id]["status"] = "Error"
        done_event.set()
        return

    updates[sra_id]["annotated_symbol"] = "[green]✓[/green]"
    updates[sra_id]["num_proteins"] = len(all_proteins) if all_proteins else None
    updates[sra_id]["status"] = "HMM Searching..."

    try:
        all_hmm_results = search_hmm_streaming(hmmer_sequences, str(hmm_file), evalue)
    except Exception:
        import traceback
        traceback.print_exc()
        updates[sra_id]["status"] = "Error"
        done_event.set()
        return

    total_hits_for_display = 0
    for hmm_name, (matching_contigs, matching_proteins) in all_hmm_results.items():
        total_hits_for_display += len(matching_contigs)

    updates[sra_id]["searched_symbol"] = "[green]✓[/green]"
    updates[sra_id]["num_hits"] = total_hits_for_display
    updates[sra_id]["status"] = "Finished"
    
    result = ProcessedSRA(
        sra_id=sra_id,
        all_contigs=all_contigs,
        all_gff=all_annotations,
        all_proteins=all_proteins,
        matching_contigs=set(),
        matching_proteins=set()
    )
    
    for hmm_name, (matching_contigs, matching_proteins) in all_hmm_results.items():
        partial_result = ProcessedSRA(
            sra_id=sra_id,
            all_contigs=all_contigs,
            all_gff=all_annotations,
            all_proteins=all_proteins,
            matching_contigs=matching_contigs,
            matching_proteins=matching_proteins
        )
        save_results(
            partial_result,
            output_dir=output_dir,
            hmm_name=hmm_name,
            save_all=save_all
        )

    del all_contigs, all_annotations, all_proteins, hmmer_sequences, result, all_hmm_results
    done_event.set()
