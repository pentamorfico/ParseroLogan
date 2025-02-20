# hmmsearch.py
from typing import List, Dict, Tuple, Set
import pyhmmer


def search_hmm_streaming(hmmer_sequences: List[pyhmmer.easel.TextSequence],
                         hmm_file: str, evalue_threshold: float) -> Dict[str, Tuple[Set[str], Set[str]]]:
    results_hmmer = {}
    matching_contigs: Set[str] = set()
    matching_proteins: Set[str] = set()

    with pyhmmer.plan7.HMMFile(hmm_file) as hmm_handle:
        for hmm in hmm_handle:
            hmm_name = hmm.name.decode()
            alphabet = hmm.alphabet
            digital_sequences = [seq.digitize(alphabet) for seq in hmmer_sequences]

            for hits in pyhmmer.hmmsearch(hmm, digital_sequences):
                for hit in hits:
                    if hit.evalue <= evalue_threshold:
                        protein_id = hit.name.decode()
                        contig_id = protein_id.rsplit("_",1)[0]
                        if protein_id not in results_hmmer:
                            results_hmmer[protein_id] = {"hmm_name":hmm_name,
                                                         "evalue":hit.evalue,
                                                         "score":hit.score,
                                                         "contig":contig_id}
                            matching_proteins.add(protein_id)
                            matching_contigs.add(contig_id)
                        else:
                            if hit.score > results_hmmer[protein_id]["score"]:
                                results_hmmer[protein_id] = {"hmm_name":hmm_name,
                                                             "evalue":hit.evalue,
                                                             "score":hit.score,
                                                             "contig":contig_id}
                        

    return results_hmmer
