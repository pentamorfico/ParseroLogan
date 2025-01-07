from typing import List, Dict, Tuple
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyrodigal_gv
import pyhmmer

from .fasta import stream_fasta_from_memory


def has_dtr(seq_record: SeqRecord, min_length: int = 30):
    substring = str(seq_record.seq).lower()[:min_length]
    pos = str(seq_record.seq).lower().rfind(substring)
    if pos < len(seq_record.seq) / 2:
        return False, 0
    substring = str(seq_record.seq).lower()[pos:]
    return str(seq_record.seq).lower()[: len(substring)] == substring, len(substring)


def fix_circle(seq_record: SeqRecord, k: int = 31):
    seq = str(seq_record.seq)
    kmers = set()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if kmer in kmers:
            new_seq = seq[: i + k - 1]
            return SeqRecord(Seq(new_seq), id=seq_record.id, description=seq_record.description)
        kmers.add(kmer)
    return seq_record


def process_contig(seq_record: SeqRecord) -> Tuple[List[Dict], List[SeqRecord], List[pyhmmer.easel.TextSequence]]:
    annotations = []
    proteins = []
    hmmer_sequences = []

    orf_finder = pyrodigal_gv.ViralGeneFinder(meta=True)
    orfs = orf_finder.find_genes(str(seq_record.seq).encode())

    for i, orf in enumerate(orfs):
        gene_id = f"{seq_record.id}_{i+1}"
        gff_entry = {
            'seqid': seq_record.id,
            'source': 'pyrodigal-gv',
            'type': 'CDS',
            'start': orf.begin + 1,
            'end': orf.end,
            'score': '.',
            'strand': '+' if orf.strand == 1 else '-',
            'phase': '0',
            'attributes': f'ID={gene_id}'
        }
        annotations.append(gff_entry)

        protein_seq = Seq(orf.translate())
        protein_record = SeqRecord(protein_seq, id=gene_id)
        proteins.append(protein_record)

        name = gene_id.encode()
        description = f"{orf.begin}..{orf.end} {orf.strand}".encode()
        sequence = orf.translate()

        hmmer_sequence = pyhmmer.easel.TextSequence(
            name=name,
            sequence=sequence,
            description=description
        )
        hmmer_sequences.append(hmmer_sequence)

    return annotations, proteins, hmmer_sequences


def process_annotations(data: bytes, min_length: int, sra_id: str, fix_circles_flag: bool = False):
    all_contigs = {}
    all_annotations = []
    all_proteins = []
    all_hmmer_sequences = []

    for contig in stream_fasta_from_memory(data, min_length, sra_id):
        dtr_check, _ = has_dtr(contig)
        if fix_circles_flag and dtr_check:
            contig = fix_circle(contig)

        annotations, proteins, hmmer_sequences = process_contig(contig)
        all_contigs[contig.id] = contig
        all_annotations.extend(annotations)
        all_proteins.extend(proteins)
        all_hmmer_sequences.extend(hmmer_sequences)

    return all_contigs, all_annotations, all_proteins, all_hmmer_sequences
