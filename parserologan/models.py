from dataclasses import dataclass
from typing import Dict, List, Set, Optional
from Bio.SeqRecord import SeqRecord

@dataclass
class ProcessedSRA:
    sra_id: str
    all_contigs: Dict[str, SeqRecord]
    all_gff: List[Dict]
    all_proteins: List[SeqRecord]
    matching_contigs: Set[str]
    matching_proteins: Set[str]
    error_message: Optional[str] = None
