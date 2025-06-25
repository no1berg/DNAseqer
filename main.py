# A light weight DNA to RNA transcription tool

# Dependencies
import create_report
import datetime
import os

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

# Store information for report (change this to a class)
sequence_info = {
    'dna_seq' : '',
    'dna_seq_len' : 0,
    'dna_gc_content' : 0.0,
    'orfs_amount' : 0,
    'orfs' : [],
    'orfs_gc' : [],
    'raw_rna_seq' : '',
    'orfs_rna' : [],
    'orfs_protein' : [],
    'protein_seq' : ''
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
# def get_dna_length(dna_sequence: str) -> int:
#     return len(dna_sequence)

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

def translate_rna_to_protein(rna_sequence: str, codon_table: dict = CODON_TABLE) -> str:
    """
    Translates an RNA sequence into a protein sequence until a stop codon is reached.

    Args:
        rna_sequence: The input RNA string (e.g., 'AUGGCUAG').
        codon_table: Optional codon table in dict structure. Default is standard.

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
        amino_acid = codon_table.get(codon, '') # Use .get() for safety against unknown codons

        # If the amino acid is an empty string (stop codon), stop translation.
        if not amino_acid:
            break
        
        protein_sequence.append(amino_acid)

    return "".join(protein_sequence)



# Define the main function
def main():
    # Set up the output directory for reports
    REPORTS_DIR =  "reports"
    os.makedirs(REPORTS_DIR, exist_ok=True)

    """Main function to run the DNA to RNA transcription tool."""
    # Prompt the user for a DNA sequence
    while True:
        try:
            dna_sequence = get_dna_sequence()
            sequence_info['dna_seq'] = dna_sequence
            print(f'Using DNA sequence: {dna_sequence}')
            break
        except ValueError as e:
            print(e)

    print('\n----- Sequence Statistics -----\n')
    # Display the length of the DNA sequence
    sequence_info['dna_seq_len'] = len(dna_sequence)
    print(f'DNA sequence length: {sequence_info["dna_seq_len"]} nucleotides\n')

    # Display the GC content of the DNA sequence
    sequence_info['dna_gc_content'] = get_gc_content(dna_sequence)
    print(f'GC content: {sequence_info["dna_gc_content"]:.2f}%\n')

    print('\n----- Open Reading Frames (ORFs) -----\n')

    # Find and display open reading frames (ORFs)
    orfs = find_all_orfs(dna_sequence)
    if orfs:
        sequence_info['orfs'] = orfs
        sequence_info['orfs_amount'] = len(orfs)
        print(f'Open Reading Frames (ORFs) found: {", ".join(orfs)}\n')
        print(f'Total number of ORFs found: {sequence_info["orfs_amount"]}\n')
        # GC content for each ORF, using a list comprehension
        orf_gc_contents = [get_gc_content(orf) for orf in orfs]
        sequence_info['orfs_gc'] = orf_gc_contents
        print(f'GC content for each ORF: {", ".join(f"{gc:.2f}%" for gc in orf_gc_contents)}\n')
    else:
        sequence_info['orfs'] = ['No ORFs present']
        sequence_info['orfs_amount'] = 0
        sequence_info['orfs_gc'] = [None]
        print('No Open Reading Frames (ORFs) found.\n')

    print('\n----- Transcription -----\n')

    # Transcribe the entire DNA sequence to RNA
    transcribed_rna = transcribe_dna_to_rna(dna_sequence)
    sequence_info['raw_rna_seq'] = transcribed_rna
    print(f'Transcribed RNA sequence: {transcribed_rna}\n')

    # Translate the RNA sequence to a protein sequence starting from the first start codon
    print("Translating raw RNA from the first start codon (AUG)...")
    start_position = transcribed_rna.find('AUG') # Find the index of the first 'AUG'

    if start_position != -1:
        # If a start codon is found, translate the sequence from that point
        protein_sequence = translate_rna_to_protein(transcribed_rna[start_position:])
        print(f"Translated protein from first ORF: {protein_sequence}\n")
    else:
        # If no start codon is found, the sequence cannot be translated
        protein_sequence = "No start codon (AUG) found."
        print("No start codon (AUG) found in the raw sequence.\n")

    sequence_info['protein_seq'] = protein_sequence

    # Transcribe the ORFs to RNA using the streamlined approach
    if orfs:
        # STREAMLINED: Reusing the transcribe_dna_to_rna function with a list comprehension
        rna_orfs = [transcribe_dna_to_rna(orf) for orf in orfs]
        print(f'Transcribed RNA from ORFs: {", ".join(rna_orfs)}\n')
        sequence_info['orfs_rna'] = rna_orfs
        # Translate the RNA sequences to protein sequences
        protein_sequences = [translate_rna_to_protein(rna) for rna in rna_orfs]
        print(f'Translated protein sequences from ORFs: {", ".join(protein_sequences)}\n')
        sequence_info['orfs_protein'] = protein_sequences
    else:
        print('No ORFs to transcribe to RNA.\n')
        sequence_info['orfs_rna'] = [None]
        sequence_info['orfs_protein'] = [None]

    #-- Report Generation --
    
    # Define a base name for the report file
    today = datetime.datetime.now().strftime("%Y-%m-%d")
    desired_base_name = os.path.join(REPORTS_DIR, f"dna_analysis_report_{today}")

    # Generate a unique base name for the report file to avoid overwriting
    unique_base_name = create_report.generate_report_filename(desired_base_name)
    print(f"Report will be saved as: {unique_base_name}.txt and {unique_base_name}.json")

    # Create the report content using the create_report module
    report_content = create_report.create_report_content(sequence_info)

    # Write the report to a file
    if create_report.write_report_to_file(report_content, unique_base_name):
        print("Report successfully created and saved.")
    else:
        print("Failed to create the report.")
    
    # Write the report to a JSON file
    if create_report.write_data_to_json(sequence_info, unique_base_name):
        print("JSON report successfully created and saved.")
    else:
        print("Failed to create the JSON report.")



# Call the main function
if __name__ == "__main__":
    main()