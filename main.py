# A light weight DNA to RNA transcription tool

# Define the standard codon table as a constant for efficiency.


CODON_TABLE = {
    # Phenylalanine (F)
    'UUU': 'F', 'UUC': 'F',
    # Leucine (L)
    'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    # Isoleucine (I)
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I',
    # Methionine (M) / Start
    'AUG': 'M',
    # Valine (V)
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    # Serine (S)
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S',
    # Proline (P)
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    # Threonine (T)
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    # Alanine (A)
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    # Tyrosine (Y)
    'UAU': 'Y', 'UAC': 'Y',
    # Histidine (H)
    'CAU': 'H', 'CAC': 'H',
    # Glutamine (Q)
    'CAA': 'Q', 'CAG': 'Q',
    # Asparagine (N)
    'AAU': 'N', 'AAC': 'N',
    # Lysine (K)
    'AAA': 'K', 'AAG': 'K',
    # Aspartic Acid (D)
    'GAU': 'D', 'GAC': 'D',
    # Glutamic Acid (E)
    'GAA': 'E', 'GAG': 'E',
    # Cysteine (C)
    'UGU': 'C', 'UGC': 'C',
    # Tryptophan (W)
    'UGG': 'W',
    # Arginine (R)
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R',
    # Glycine (G)
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    # Stop Codons (translates to an empty string to terminate translation)
    'UAA': '', 'UAG': '', 'UGA': ''
}


# Request the user for a DNA sequence
def get_dna_sequence() -> str:
    dna_sequence = input("Please enter a DNA sequence (A, T, C, G): ").upper()
    if not all(base in "ATCG" for base in dna_sequence):
        raise ValueError("Invalid DNA sequence. Please use only A, T, C, and G.")
    return dna_sequence.upper()

# Generate the reverse complement of a DNA sequence
def reverse_complement(dna_sequence: str) -> str:
    """Generates the reverse complement of a DNA sequence."""
    complement_map = str.maketrans('ATCG', 'TAGC')
    return dna_sequence.upper().translate(complement_map)[::-1]

# Find the length of the DNA sequence
def get_dna_length(dna_sequence: str) -> int:
    return len(dna_sequence)

# Find the GC content of the DNA sequence
def get_gc_content(dna_sequence: str) -> float:
    if not dna_sequence:
        return 0.0
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    # Calculate GC content as a percentage
    gc_content = (gc_count / len(dna_sequence)) * 100
    return gc_content

# Transcribe a DNA sequence to an RNA sequence
def transcribe_dna_to_rna(dna_sequence: str) -> str:
    """Transcribes a single DNA string to an RNA string."""
    return dna_sequence.replace('T', 'U')

# Find all open reading frames (ORFs)
def find_all_orfs(dna_sequence: str, min_orf_length: int = 50) -> list:
    """Find all open reading frames (ORFs) in the six reading frames of a DNA sequence.
    Args:
        dna_sequence: The input DNA string (e.g., 'ATGGCTAG').
        min_orf_length: Minimum length of ORFs to consider (default is 50)."""
    all_orfs = []
    strands = {
        'forward': dna_sequence,
        'reverse': reverse_complement(dna_sequence)
    }

    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']

    for strand_name, sequence in strands.items():
        for frame in range(3): # Corresponds to starting at positions 0, 1, and 2
            # Find all the start and stop codon positions in the current reading frame
            starts = [i for i in range(frame, len(sequence), 3) if sequence[i:i+3] == start_codon]
            stops = [i for i in range(frame, len(sequence), 3) if sequence[i:i+3] in stop_codons]

            # For each start codon, find the first stop codon that follows
            for start_pos in starts:
                for stop_pos in stops:
                    if stop_pos > start_pos:
                        orf = sequence[start_pos:stop_pos+3]
                        if len(orf) >= min_orf_length:
                            all_orfs.append(orf)
                        break  # Stop after the first valid stop codon for this start codon
    return all_orfs

def translate_rna_to_protein(rna_sequence: str) -> str:
    """
    Translates an RNA sequence into a protein sequence until a stop codon is reached.

    Args:
        rna_sequence: The input RNA string (e.g., 'AUGGCUAG').

    Returns:
        The resulting protein sequence as a string (e.g., 'MA').
    """
    protein_sequence = []
    # Iterate through the RNA sequence in steps of 3 (codon length)
    for i in range(0, len(rna_sequence), 3):
        # Ensure we have a full codon to translate
        if i + 3 > len(rna_sequence):
            break
            
        codon = rna_sequence[i:i+3]
        amino_acid = CODON_TABLE.get(codon, '') # Use .get() for safety against unknown codons

        # If the amino acid is an empty string (stop codon), stop translation.
        if not amino_acid:
            break
        
        protein_sequence.append(amino_acid)

    return "".join(protein_sequence)




# Define the main function
def main():
    """Main function to run the DNA to RNA transcription tool."""
    # Prompt the user for a DNA sequence
    while True:
        try:
            dna_sequence = get_dna_sequence()
            print(f'Using DNA sequence: {dna_sequence}')
            break
        except ValueError as e:
            print(e)

    print('\n----- Sequence Statistics -----\n')
    # Display the length of the DNA sequence
    print(f'DNA sequence length: {get_dna_length(dna_sequence)} nucleotides\n')

    # Display the GC content of the DNA sequence
    print(f'GC content: {get_gc_content(dna_sequence):.2f}%\n')

    print('\n----- Open Reading Frames (ORFs) -----\n')

    # Find and display open reading frames (ORFs)
    orfs = find_all_orfs(dna_sequence)
    if orfs:
        print(f'Open Reading Frames (ORFs) found: {", ".join(orfs)}\n')
        print(f'Total number of ORFs found: {len(orfs)}\n')
        # GC content for each ORF, using a list comprehension
        orf_gc_contents = [get_gc_content(orf) for orf in orfs]
        print(f'GC content for each ORF: {", ".join(f"{gc:.2f}%" for gc in orf_gc_contents)}\n')
    else:
        print('No Open Reading Frames (ORFs) found.\n')

    print('\n----- Transcription -----\n')

    # Transcribe the entire DNA sequence to RNA
    transcribed_rna = transcribe_dna_to_rna(dna_sequence)
    print(f'Transcribed RNA sequence: {transcribed_rna}\n')

    # Transcribe the ORFs to RNA using the streamlined approach
    if orfs:
        # STREAMLINED: Reusing the transcribe_dna_to_rna function with a list comprehension
        rna_orfs = [transcribe_dna_to_rna(orf) for orf in orfs]
        print(f'Transcribed RNA from ORFs: {", ".join(rna_orfs)}\n')
        # Translate the RNA sequences to protein sequences
        protein_sequences = [translate_rna_to_protein(rna) for rna in rna_orfs]
        print(f'Translated protein sequences from ORFs: {", ".join(protein_sequences)}\n')
    else:
        print('No ORFs to transcribe to RNA.\n')

# Call the main function
if __name__ == "__main__":
    main()