# Rosalind Solutions

A small, dependency-free Python library solving ten problems from the [Rosalind Bioinformatics Stronghold](https://rosalind.info/problems/list-view/). All solutions are in **`rosalind_functions.py`**. Functions return values instead of printing them.

## Problems and functions

| Problem | Description | Function | Return value |
|---|---|---|---|
| [DNA](https://rosalind.info/problems/dna/) | Count the four DNA nucleotides. | `count_dna_nucleotides(dna)` | Tuple `(A, C, G, T)` |
| [RNA](https://rosalind.info/problems/rna/) | Transcribe DNA by changing T to U. | `transcribe_dna(dna)` | RNA string |
| [REVC](https://rosalind.info/problems/revc/) | Find the reverse complement of DNA. | `reverse_complement(dna)` | DNA string |
| [GC](https://rosalind.info/problems/gc/) | Calculate GC percentage for a DNA sequence. | `gc_content(dna)` | Percentage as `float` |
| [HAMM](https://rosalind.info/problems/hamm/) | Count nucleotide differences at aligned positions. | `hamming_distance(first, second)` | Integer |
| [PROT](https://rosalind.info/problems/prot/) | Translate an RNA sequence into a protein, stopping at a stop codon. | `translate_rna(rna)` | Protein string |
| [SUBS](https://rosalind.info/problems/subs/) | Find all occurrences of a DNA motif, including overlaps. | `find_motif(dna, motif)` | List of 1-based positions |
| [PRTM](https://rosalind.info/problems/prtm/) | Sum monoisotopic masses of amino acid residues. | `protein_mass(protein)` | Mass in Da as `float` |
| [REVP](https://rosalind.info/problems/revp/) | Find reverse-complement palindromes of length 4–12. | `reverse_palindromes(dna)` | List of `(1-based position, length)` tuples |
| [TRAN](https://rosalind.info/problems/tran/) | Compute transition/transversion ratio between aligned DNA strings. | `transition_transversion_ratio(first, second)` | Ratio as `float` |

## Usage

Requires Python 3; no third-party libraries.

```python
from rosalind_functions import count_dna_nucleotides, find_motif, translate_rna

print(count_dna_nucleotides('AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC'))
print(find_motif('GATATATGCATATACTT', 'ATAT'))
print(translate_rna('AUGGCCUAA'))
```

The functions themselves never print. The `print` calls above only demonstrate their return values.

## Notes

- DNA and RNA input strings should be uppercase and contain the standard nucleotide letters. `GC` operates on a single DNA string; the Rosalind GC problem supplies multiple FASTA records, so call it on each record's sequence and select the largest percentage separately.
- `HAMM` and `TRAN` require equal-length aligned strings. `TRAN` raises `ZeroDivisionError` when there are no transversions.
- `SUBS` and `REVP` use **1-based positions**, as Rosalind does.
- `PRTM` uses Rosalind's **residue masses** (no extra water mass).
- `REVP` takes the DNA string, not a FASTA header. Read and concatenate FASTA sequence lines before passing it in.
- Rosalind expects space-separated output for some problems. The functions return structured Python values instead; format them when submitting to Rosalind.

## Run tests

```bash
python3 test_rosalind.py
```

## Git history

The project was developed in five stages: basic DNA operations, translation/motifs, remaining bioinformatics functions, documentation, and tests.
