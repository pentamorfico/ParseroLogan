from pathlib import Path
from .models import ProcessedSRA

def save_results(processed_sra: ProcessedSRA, output_dir: Path, hmm_name: str = None, save_all: bool = False):
    if not (processed_sra.matching_contigs or save_all):
        print(f"No matching contigs and save_all is False. No files will be saved for {processed_sra.sra_id}.")
        return

    if hmm_name is None:
        sra_output_dir = output_dir / processed_sra.sra_id
    else:
        sra_output_dir = output_dir / hmm_name / processed_sra.sra_id

    sra_output_dir.mkdir(parents=True, exist_ok=True)

    if save_all:
        with open(sra_output_dir / f"{processed_sra.sra_id}.fna", 'w') as f:
            for contig in processed_sra.all_contigs.values():
                f.write(f">{contig.description}\n{contig.seq}\n")

        with open(sra_output_dir / f"{processed_sra.sra_id}.faa", 'w') as f:
            for protein in processed_sra.all_proteins:
                f.write(f">{protein.id}\n{protein.seq}\n")

        with open(sra_output_dir / f"{processed_sra.sra_id}.gff", 'w') as f:
            f.write("##gff-version 3\n")
            for ann in processed_sra.all_gff:
                gff_line = '\t'.join(str(ann[key]) for key in
                                     ['seqid', 'source', 'type', 'start', 'end',
                                      'score', 'strand', 'phase', 'attributes'])
                f.write(gff_line + "\n")

        if processed_sra.matching_proteins:
            with open(sra_output_dir / f"{processed_sra.sra_id}_protein_hits.txt", 'w') as f:
                for pid in processed_sra.matching_proteins:
                    f.write(pid + "\n")
    else:
        with open(sra_output_dir / f"{processed_sra.sra_id}.fna", 'w') as f:
            for contig_id in processed_sra.matching_contigs:
                if contig_id in processed_sra.all_contigs:
                    contig = processed_sra.all_contigs[contig_id]
                    f.write(f">{contig.description}\n{contig.seq}\n")

        with open(sra_output_dir / f"{processed_sra.sra_id}.faa", 'w') as f:
            for protein in processed_sra.all_proteins:
                cid = protein.id.rsplit('_', 1)[0]
                if cid in processed_sra.matching_contigs:
                    f.write(f">{protein.id}\n{protein.seq}\n")

        with open(sra_output_dir / f"{processed_sra.sra_id}.gff", 'w') as f:
            f.write("##gff-version 3\n")
            for ann in processed_sra.all_gff:
                if ann['seqid'] in processed_sra.matching_contigs:
                    gff_line = '\t'.join(str(ann[key]) for key in
                                         ['seqid', 'source', 'type', 'start', 'end',
                                          'score', 'strand', 'phase', 'attributes'])
                    f.write(gff_line + "\n")

        if processed_sra.matching_proteins:
            with open(sra_output_dir / f"{processed_sra.sra_id}_protein_hits.txt", 'w') as f:
                for pid in processed_sra.matching_proteins:
                    f.write(pid + "\n")
