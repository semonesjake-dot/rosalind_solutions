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
