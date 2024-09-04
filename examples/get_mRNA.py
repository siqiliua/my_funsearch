import itertools
import numpy as np
import subprocess

# Example of codon table for human genome
CODON_TABLE = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'C': ['CUU', 'UUC', 'UUA', 'UUG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['UUU', 'UUC'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'K': ['AAA', 'AAG'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CCA', 'CCG'],
    'M': ['AUG'],
    'N': ['AAU', 'AAC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAG', 'CAA'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    '*': ['UAA', 'UAG', 'UGA']  # Stop codon
}


def run_rnafold(sequence):
    try:
        process = subprocess.Popen(
            ['RNAfold'],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        stdout, stderr = process.communicate(sequence.encode())

        if process.returncode != 0:
            print("Error:", stderr.decode())
            return None

        output = stdout.decode().strip()
        lines = output.split('\n')
        structure, free_energy = lines[1].split(' (')
        free_energy = float(free_energy.rstrip(')'))

        return structure, free_energy

    except FileNotFoundError:
        print("RNAfold is not installed or not found in PATH.")
        return None


def translate_codon(codon):
    """Translate a codon into its corresponding amino acid."""
    for amino_acid, codons in CODON_TABLE.items():
        if codon in codons:
            return amino_acid
    return None

@funsearch.run
def evaluate(protein_sequence: str) -> str:
    """Evaluates the protein sequence and returns the corresponding base sequence."""
    base_sequence = solve(protein_sequence)
    structure, energy = run_rnafold(base_sequence)
    print(f"Structure: {structure}\nFree Energy: {energy} kcal/mol")
    return energy


# def solve(protein_sequence: str) -> str:
#     """Constructs a valid RNA sequence from a protein sequence."""
#     rna_sequence = ''
#     for amino_acid in protein_sequence:
#         codons = CODON_TABLE.get(amino_acid)
#         if codons:
#             rna_sequence += priority(codons)  # Selecting the first codon for simplicity
#         else:
#             raise ValueError(f"Unknown amino acid: {amino_acid}")
#     return rna_sequence

def solve(protein_sequence: str) -> str:
    """Constructs a valid RNA sequence from a protein sequence."""
    rna_sequence = ''
    for amino_acid in protein_sequence:
        codons = CODON_TABLE.get(amino_acid)
        if codons:
            # Using priority to choose the codon with the highest priority score
            priorities = np.array([priority(codon) for codon in codons])
            best_codon_index = np.argmax(priorities)  # Get the index of the highest priority codon
            rna_sequence += codons[best_codon_index]  # Select the codon with highest priority
        else:
            raise ValueError(f"Unknown amino acid: {amino_acid}")
    return rna_sequence


@funsearch.evolve
def priority(codon: str) -> float:
    """Define a priority scoring for codons. Custom logic should be implemented here."""
    # Example: Higher scores for codons that are more frequently used
    # This is a simple example. You can modify it based on your criteria.

    return 0.0


# Example usage
# protein_seq = "ACGT"
# rna_seq = evaluate(protein_seq)
# print(f"Corresponding RNA Sequence: {rna_seq}")
