from io import BytesIO, TextIOWrapper
from typing import Iterator
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def stream_fasta_from_memory(data: bytes, min_length: int, sra_id: str = "") -> Iterator[SeqRecord]:
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


def count_fasta_contigs(data: bytes) -> int:
    text_stream = TextIOWrapper(BytesIO(data))
    count = 0
    for line in text_stream:
        if line.startswith('>'):
            count += 1
    return count
