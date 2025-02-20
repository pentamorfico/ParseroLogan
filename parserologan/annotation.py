# annotation.py
from typing import List, Dict, Tuple, Iterator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyrodigal_gv
import pyhmmer
from io import BytesIO, TextIOWrapper

def filter_fasta(data: bytes, min_length: int, sra_id: str = "") -> Iterator[SeqRecord]:
    text_stream = TextIOWrapper(BytesIO(data))
    current_header = ""
    current_seq = []

    def yield_contig():
        seq = ''.join(current_seq)
        if len(seq) >= min_length:
            seq_id = current_header[1:].split()[0]
            description = current_header[1:]
            if seq_id.startswith("_") and sra_id:
                seq_id = f"{sra_id}{seq_id}"
                description = f"{sra_id}{current_header[1:]}"
            return SeqRecord(Seq(seq), id=seq_id, description=description)
        return None

    for line in text_stream:
        line = line.strip()
        if line.startswith('>'):
            if current_seq:
                contig = yield_contig()
                if contig:
                    yield contig
            current_header = line
            current_seq = []
        else:
            current_seq.append(line)
    if current_seq:
        contig = yield_contig()
        if contig:
            yield contig

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

def process_annotations(data: bytes, min_length: int, sra_id: str):
    all_contigs = {}
    all_annotations = []
    all_proteins = []
    all_hmmer_sequences = []

    for contig in filter_fasta(data, min_length, sra_id):
        annotations, proteins, hmmer_sequences = process_contig(contig)
        all_contigs[contig.id] = contig
        all_annotations.extend(annotations)
        all_proteins.extend(proteins)
        all_hmmer_sequences.extend(hmmer_sequences)

    return all_contigs, all_annotations, all_proteins, all_hmmer_sequences
