# download.py
import subprocess
from pathlib import Path

from .annotation import process_annotations
from .hmm_search import search_hmm_streaming
from io import BytesIO, TextIOWrapper

def count_fasta_contigs(data: bytes) -> int:
    
    text_stream = TextIOWrapper(BytesIO(data))
    count = 0
    for line in text_stream:
        if line.startswith('>'):
            count += 1
    return count
# Create a global variable to hold the manager.
global_manager = None

def set_manager(mgr):
    """Set the global manager. This must be called once in the main process before spawning workers."""
    global global_manager
    global_manager = mgr


def process_sample(sra_id: str, updates: dict, done_event, hmm_file: Path, min_length: int, evalue: float, output_dir: Path, save_all: bool = False):
    # Use the global manager to create a proxy dict for this SRA.
    updates[sra_id] = global_manager.dict({
        "name": sra_id,
        "status": "Downloading...",
        "downloaded_symbol": "[yellow]?[/yellow]",
        "annotated_symbol": "",
        "searched_symbol": "",
        "num_contigs": None,
        "num_proteins": None,
        "num_hits": None
    })
    bucket_key = f"s3://logan-pub/c/{sra_id}/{sra_id}.contigs.fa.zst"
    cmd = f"aws s3 cp {bucket_key} - --no-sign-request --no-progress | seqkit seq -m {min_length}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    try:
        data = process.stdout.read()
    except Exception:
        print("failed at download")
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
    num_contigs = count_fasta_contigs(data)
    updates[sra_id]["num_contigs"] = num_contigs if num_contigs > 0 else None
    updates[sra_id]["status"] = "Annotating..."

    try:
        all_contigs, all_annotations, all_proteins, hmmer_digital_seqs = process_annotations(
            data, min_length, sra_id
        )
    except Exception:
        import traceback
        print("failed pyrodigal")
        traceback.print_exc()
        updates[sra_id]["status"] = "Error"
        done_event.set()
        return

    updates[sra_id]["annotated_symbol"] = "[green]✓[/green]"
    updates[sra_id]["num_proteins"] = len(all_proteins) if all_proteins else None
    updates[sra_id]["status"] = "HMM Searching..."
    
    try:
        all_hmm_results = search_hmm_streaming(hmmer_digital_seqs, str(hmm_file), evalue)
    except Exception:
        import traceback
        print("failed hmmsearch")
        traceback.print_exc()
        updates[sra_id]["status"] = "Error"
        done_event.set()
        return

    total_hits_for_display = 0
    total_hits_for_display += len(all_hmm_results)

    updates[sra_id]["searched_symbol"] = "[green]✓[/green]"
    updates[sra_id]["num_hits"] = total_hits_for_display
    updates[sra_id]["status"] = "Finished"
    
    matching_contigs = [all_hmm_results[n]["contig"] for n in all_hmm_results]
    matching_proteins = [n for n in all_hmm_results]
    
    if not save_all:
        if len(matching_proteins)>0:
            saved_contigs = {contig_id: contig for contig_id, contig in all_contigs.items() if contig_id in matching_contigs}
            saved_annots = [annot for annot in all_annotations if annot["seqid"] in matching_contigs]
            saved_proteins = [n for n in all_proteins if n.id.rsplit('_', 1)[0] in matching_contigs]
        else:
            print(f"No matching contigs and save_all is False. No files will be saved for {sra_id}.")
            return
    else:
        saved_contigs, saved_annots, saved_proteins = all_contigs, all_annotations, all_proteins
        
    sra_output_dir = output_dir / sra_id
    sra_output_dir.mkdir(parents=True, exist_ok=True)
    

        
    with open(sra_output_dir / f"{sra_id}.fna", 'w') as fna, \
        open(sra_output_dir / f"{sra_id}.faa", 'w') as faa, \
        open(sra_output_dir / f"{sra_id}_prot_hits.txt", 'w') as hitfile, \
        open(sra_output_dir / f"{sra_id}.gff", 'w') as gff:
            
            print(f"Saving results for {sra_id}.")

            
            for contig in saved_contigs.values():
                fna.write(f">{contig.description}\n{contig.seq}\n")
            
            for ann in saved_annots:
                gff_line = '\t'.join(str(ann[key]) for key in
                        ['seqid', 'source', 'type', 'start', 'end',
                        'score', 'strand', 'phase', 'attributes'])
                gff.write(gff_line + "\n")
                
            for protein in saved_proteins:
                faa.write(f">{protein.id}\n{protein.seq}\n")
                
            for hmm_hit in all_hmm_results:
                hmm_hit_dicc = all_hmm_results[hmm_hit]
                #write tab-separated with the protein id, hmm_name, evalue, score, contig
                hitfile.write(f"{hmm_hit}\t{hmm_hit_dicc['hmm_name']}\t{hmm_hit_dicc['evalue']}\t{hmm_hit_dicc['score']}\t{hmm_hit_dicc['contig']}\n")
                

    
    done_event.set()
