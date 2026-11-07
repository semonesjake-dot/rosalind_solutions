"""Simple tests using Rosalind sample inputs; no external test framework."""
from rosalind_functions import (
    count_dna_nucleotides, transcribe_dna, reverse_complement, gc_content,
    hamming_distance, translate_rna, find_motif, protein_mass,
    reverse_palindromes, transition_transversion_ratio,
)


def run_tests():
    assert count_dna_nucleotides('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC') == (20, 12, 17, 21)
    assert transcribe_dna('GATGGAACTTGACTACGTAAATT') == 'GAUGGAACUUGACUACGUAAAUU'
    assert reverse_complement('AAAACCCGGT') == 'ACCGGGTTTT'
    assert gc_content('AGCT') == 50.0
    assert hamming_distance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT') == 7
    assert translate_rna('AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG') == 'MAIVMGR'
    assert find_motif('GATATATGCATATACTT', 'ATAT') == [2, 4, 10]
    assert abs(protein_mass('SKADYEK') - 821.392) < 0.001
    assert reverse_palindromes('TCAATGCATGCGGGTCTATATGCAT') == [(4, 6), (5, 4), (6, 6), (7, 4), (17, 4), (18, 4), (20, 6), (21, 4)]
    assert transition_transversion_ratio('AG', 'GC') == 1.0
    print('All 10 Rosalind function tests passed.')


if __name__ == '__main__':
    run_tests()
