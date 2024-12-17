from dataclasses import dataclass
from typing import Dict, List, Set, Optional
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def count_fasta_contigs(data: bytes) -> int:
    from io import BytesIO, TextIOWrapper
    text_stream = TextIOWrapper(BytesIO(data))
    count = 0
    for line in text_stream:
        if line.startswith('>'):
            count += 1
    return count

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

def display_value(val):
    return str(val) if val is not None else ""

@dataclass
class ProcessedSRA:
    sra_id: str
    all_contigs: Dict[str, SeqRecord]
    all_gff: List[Dict]
    all_proteins: List[SeqRecord]
    matching_contigs: Set[str]
    matching_proteins: Set[str]
    error_message: Optional[str] = None
