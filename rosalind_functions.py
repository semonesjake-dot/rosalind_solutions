"""Solutions to ten Rosalind Bioinformatics Stronghold problems.

Each function returns its answer without printing or external dependencies.
Positions returned by motif searches use Rosalind's 1-based indexing.
"""


def count_dna_nucleotides(dna):
    """DNA: Return counts of A, C, G, and T in that order as a tuple."""
    return tuple(dna.count(base) for base in 'ACGT')


def transcribe_dna(dna):
    """RNA: Return the RNA transcript of a coding-strand DNA sequence."""
    return dna.replace('T', 'U')


def reverse_complement(dna):
    """REVC: Return the reverse complement of a DNA sequence."""
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complements[base] for base in reversed(dna))


def gc_content(dna):
    """GC: Return the GC percentage of a nonempty DNA sequence."""
    return 100.0 * (dna.count('G') + dna.count('C')) / len(dna)


def hamming_distance(first, second):
    """HAMM: Count differing positions in two equal-length DNA strings."""
    if len(first) != len(second):
        raise ValueError('Sequences must have equal lengths')
    return sum(a != b for a, b in zip(first, second))


# Standard RNA genetic code. '*' represents a stop codon.
RNA_CODON_TABLE = {
    'UUU':'F','UUC':'F','UUA':'L','UUG':'L',
    'UCU':'S','UCC':'S','UCA':'S','UCG':'S',
    'UAU':'Y','UAC':'Y','UAA':'*','UAG':'*',
    'UGU':'C','UGC':'C','UGA':'*','UGG':'W',
    'CUU':'L','CUC':'L','CUA':'L','CUG':'L',
    'CCU':'P','CCC':'P','CCA':'P','CCG':'P',
    'CAU':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'CGU':'R','CGC':'R','CGA':'R','CGG':'R',
    'AUU':'I','AUC':'I','AUA':'I','AUG':'M',
    'ACU':'T','ACC':'T','ACA':'T','ACG':'T',
    'AAU':'N','AAC':'N','AAA':'K','AAG':'K',
    'AGU':'S','AGC':'S','AGA':'R','AGG':'R',
    'GUU':'V','GUC':'V','GUA':'V','GUG':'V',
    'GCU':'A','GCC':'A','GCA':'A','GCG':'A',
    'GAU':'D','GAC':'D','GAA':'E','GAG':'E',
    'GGU':'G','GGC':'G','GGA':'G','GGG':'G',
}


def translate_rna(rna):
    """PROT: Translate RNA until the first stop codon (not included)."""
    protein = ''
    for i in range(0, len(rna) - 2, 3):
        amino_acid = RNA_CODON_TABLE[rna[i:i + 3]]
        if amino_acid == '*':
            break
        protein += amino_acid
    return protein


def find_motif(dna, motif):
    """SUBS: Return all 1-based starting positions, including overlaps."""
    if not motif:
        raise ValueError('Motif must not be empty')
    return [i + 1 for i in range(len(dna) - len(motif) + 1)
            if dna[i:i + len(motif)] == motif]
