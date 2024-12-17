from typing import List, Set, Tuple
import pyhmmer

def search_hmm_streaming(hmmer_sequences: List[pyhmmer.easel.TextSequence], hmm_file: str, evalue_threshold: float) -> Tuple[Set[str], Set[str]]:
    matching_contigs = set()
    matching_proteins = set()

    with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
        hmm = next(hmm_handle)
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

    return matching_contigs, matching_proteins
