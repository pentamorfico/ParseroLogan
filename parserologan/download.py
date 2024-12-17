import json
import subprocess
import threading
from pathlib import Path

from .models import ProcessedSRA
from .annotation import process_annotations
from .hmm_search import search_hmm_streaming
from .output import save_results
from .fasta import count_fasta_contigs

def get_s3_object_size(sra_id: str) -> str:
    bucket = "logan-pub"
    key = f"c/{sra_id}/{sra_id}.contigs.fa.zst"
    cmd = [
        "aws",
        "s3api",
        "head-object",
        "--bucket",
        bucket,
        "--key",
        key,
        "--no-sign-request"
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        info = json.loads(result.stdout)
        size_bytes = info["ContentLength"]
        size_mb = size_bytes / (1024 * 1024)
        return f"{size_mb:.2f} MB"
    except Exception:
        return "?"

def download_sra_to_variable(sra_id: str, updates: dict, done_event: threading.Event, hmm_file: Path, min_length: int, evalue: float, output_dir: Path, fix_circles_flag: bool = False):
    total_size = get_s3_object_size(sra_id)
    updates[sra_id] = {
        "name": sra_id,
        "total_size": total_size,
        "status": "Downloading...",
        "downloaded_symbol": "[yellow]?[/yellow]",
        "annotated_symbol": "",
        "searched_symbol": "",
        "data": b"",
        "num_contigs": None,
        "num_proteins": None,
        "num_hits": None
    }

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
            updates[sra_id]["data"], min_length, fix_circles_flag=fix_circles_flag
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
        matching_contigs, matching_proteins = search_hmm_streaming(hmmer_sequences, str(hmm_file), evalue)
    except Exception:
        import traceback
        traceback.print_exc()
        updates[sra_id]["status"] = "Error"
        done_event.set()
        return

    num_hits = len(matching_contigs) if matching_contigs else 0

    updates[sra_id]["searched_symbol"] = "[green]✓[/green]"
    updates[sra_id]["num_hits"] = num_hits
    updates[sra_id]["status"] = "Finished"

    result = ProcessedSRA(
        sra_id=sra_id,
        all_contigs=all_contigs,
        all_gff=all_annotations,
        all_proteins=all_proteins,
        matching_contigs=matching_contigs,
        matching_proteins=matching_proteins
    )
    save_results(result, output_dir=output_dir, save_all=False)

    del all_contigs, all_annotations, all_proteins, hmmer_sequences, matching_contigs, matching_proteins, result

    done_event.set()
