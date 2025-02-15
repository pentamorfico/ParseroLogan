# hmmsearch.py
from typing import List, Dict, Tuple, Set
import pyhmmer


def search_hmm_streaming(hmmer_sequences: List[pyhmmer.easel.TextSequence], hmm_file: str, evalue_threshold: float) -> Dict[str, Tuple[Set[str], Set[str]]]:
    results_per_hmm: Dict[str, Tuple[Set[str], Set[str]]] = {}

    with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
        for hmm in hmm_handle:
            hmm_name = hmm.name.decode()
            matching_contigs: Set[str] = set()
            matching_proteins: Set[str] = set()

            alphabet = hmm.alphabet
            digital_sequences = [seq.digitize(alphabet) for seq in hmmer_sequences]

            for hits in pyhmmer.hmmsearch(hmm, digital_sequences):
                for hit in hits:
                    if hit.evalue <= evalue_threshold:
                        protein_id = hit.name.decode()
                        matching_proteins.add(protein_id)
                        hit_parts = protein_id.split('_')
                        if len(hit_parts) >= 2:
                            contig_id = f"{hit_parts[0]}_{hit_parts[1]}"
                            matching_contigs.add(contig_id)

            results_per_hmm[hmm_name] = (matching_contigs, matching_proteins)

    return results_per_hmm
