# A light weight DNA to RNA transcription tool

# Dependencies
import json
import create_report
import datetime
import os
from typing import TypeVar, Dict, Any

# TypeVar to create a generic type
_T = TypeVar("_T")

# The standard codon table
CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L',
    'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UCU': 'S', 'UCC': 'S',
    'UCA': 'S', 'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'CCU': 'P', 'CCC': 'P',
    'CCA': 'P', 'CCG': 'P', 'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
    'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R',
    'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'UAA': '', 'UAG': '', 'UGA': ''
}

class SequenceInfo:
    """
    Holds and processes DNA, RNA, and protein sequence information.
    This class is now serializable to and from a dictionary.
    """
    def __init__(self, dna_seq: str = ''):
        self.dna_seq: str = dna_seq
        self.header: str = ''
        self.dna_seq_len: int = len(dna_seq)
        self.dna_gc_content: float = self.calculate_gc_content()
        self.raw_rna_seq: str = self.transcribe_dna_to_rna(self.dna_seq)
        self.orfs: list[str] = self.find_all_orfs()
        self.orfs_amount: int = len(self.orfs)
        self.orfs_gc: list[float] = [self.calculate_gc_content(orf) for orf in self.orfs]
        self.orfs_rna: list[str] = [self.transcribe_dna_to_rna(orf) for orf in self.orfs]
        self.orfs_protein: list[str] = [self.translate_rna_to_protein(rna) for rna in self.orfs_rna]
        self.protein_seq: str = self.translate_rna_to_protein(self.raw_rna_seq)

    def to_dict(self) -> Dict[str, Any]:
        """Serializes the object attributes to a dictionary."""
        return self.__dict__

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "SequenceInfo":
        """Creates a SequenceInfo instance from a dictionary."""
        instance = cls(data.get('dna_seq', ''))
        for key, value in data.items():
            setattr(instance, key, value)
        return instance

    def reverse_complement(self) -> str:
        """Generates the reverse complement of the instance's DNA sequence."""
        if not self.dna_seq:
            return ""
        complement_map = str.maketrans('ATCG', 'TAGC')
        return self.dna_seq.upper().translate(complement_map)[::-1]

    def calculate_gc_content(self, sequence: str | None = None) -> float:
        """Calculates GC content of a given sequence or the instance's DNA sequence."""
        seq = self.dna_seq if sequence is None else sequence
        if not seq:
            return 0.0
        gc_count = seq.count('G') + seq.count('C')
        return (gc_count / len(seq)) * 100

    def transcribe_dna_to_rna(self, sequence: str | None = None) -> str:
        """Transcribes a given DNA sequence to RNA."""
        seq_to_transcribe = self.dna_seq if sequence is None else sequence
        if not seq_to_transcribe:
            return ''
        return seq_to_transcribe.replace('T', 'U')

    def find_all_orfs(self, min_orf_length: int = 50) -> list[str]:
        """Finds all open reading frames (ORFs) in the six reading frames."""
        all_orfs = []
        strands = {'forward': self.dna_seq, 'reverse': self.reverse_complement()}
        start_codon = 'ATG'
        stop_codons = ['TAA', 'TAG', 'TGA']

        for sequence in strands.values():
            for frame in range(3):
                starts = [i for i in range(frame, len(sequence), 3) if sequence[i:i+3] == start_codon]
                stops = [i for i in range(frame, len(sequence), 3) if sequence[i:i+3] in stop_codons]
                for start_pos in starts:
                    for stop_pos in stops:
                        if stop_pos > start_pos:
                            orf = sequence[start_pos:stop_pos + 3]
                            if len(orf) >= min_orf_length:
                                all_orfs.append(orf)
                            break
        return all_orfs

    def translate_rna_to_protein(self, rna_sequence: str, codon_table: dict | None = None) -> str:
        """
        Translates an RNA sequence into a protein sequence, starting from the first
        start codon (AUG) until a stop codon is reached.
        """
        if codon_table is None:
            codon_table = CODON_TABLE

        # Find the first start codon 'AUG'
        start_codon_pos = rna_sequence.find('AUG')
        if start_codon_pos == -1:
            return ""  # No start codon found, so no protein

        protein_sequence = []
        # Start translation from the found start codon
        for i in range(start_codon_pos, len(rna_sequence) - 2, 3):
            codon = rna_sequence[i:i+3]
            amino_acid = codon_table.get(codon, '')
            
            # Stop translation if a stop codon is found
            if not amino_acid:
                break
            
            protein_sequence.append(amino_acid)
            
        return "".join(protein_sequence)

def get_dna_sequence() -> str:
    """Requests and validates a DNA sequence from the user."""
    while True:
        dna_sequence = input("Please enter a DNA sequence (A, T, C, G): ").upper()
        if all(base in "ATCG" for base in dna_sequence):
            return dna_sequence
        print("Invalid DNA sequence. Please use only A, T, C, and G.")

# Add function to import several sequences and header from a fasta file
def import_sequence_from_fasta(file_path: str) -> list[tuple[str, str]]:
    """
    Imports DNA sequences and their headers from a FASTA file.
    and removes any whitespace or newline characters.

    Args:
        file_path (str): The path to the FASTA file.

    Returns:
        tuple[str, str]: A tuple containing the header and the DNA sequence.
    """
    with open(file_path, 'r') as fasta_file:
        fasta_content = fasta_file.read()
        # Extract the header (first line) and sequence (remove whitespace)
        headers = []
        sequences = []
        for line in fasta_content.splitlines():
            line = line.strip()
            if line.startswith(">"):
                headers.append(line[1:])  # Remove '>' from header
            else:
                sequences.append(line)
        # Combine headers and sequences into a list of tuples
        return list(zip(headers, sequences))

def main():
    """Main function to run the DNA analysis tool."""
    REPORTS_DIR = "reports"
    os.makedirs(REPORTS_DIR, exist_ok=True)

    sequences_to_process = []

    try:
        # --- PHASE 1: Data Collection ---
        input_method = input("Choose input method (1: Manual, 2: Fasta file): ").strip()
        
        if input_method == '2':
            # FASTA file input: process multiple sequences
            fasta_file_path = input("Enter the path to the FASTA file: ").strip()
            if not os.path.isfile(fasta_file_path):
                raise FileNotFoundError(f"The specified file does not exist: {fasta_file_path}")
            
            # import_sequence_from_fasta should return a list of (header, sequence) tuples
            imported_sequences = import_sequence_from_fasta(fasta_file_path)
            if not imported_sequences:
                print("FASTA file is empty or could not be parsed.")
                return

            for header, dna_input in imported_sequences:
                sequence_info = SequenceInfo(dna_input)
                sequence_info.header = header
                sequences_to_process.append(sequence_info)
            
        elif input_method == '1':
            # Manual input: process a single sequence
            dna_input = get_dna_sequence()
            sequence_info = SequenceInfo(dna_input)
            sequence_info.header = "Manual_Input"
            sequences_to_process.append(sequence_info)
            
        else:
            raise ValueError("Invalid input method selected. Please choose '1' or '2'.")

        # --- PHASE 2: Analysis and Reporting for each sequence ---
        if not sequences_to_process:
            print("No sequences to analyze.")
            return

        print(f"\nFound {len(sequences_to_process)} sequence(s) to analyze.")

        for i, sequence_info in enumerate(sequences_to_process, start=1):
            print(f"\n\n--- Analyzing Sequence {i} of {len(sequences_to_process)}: '{sequence_info.header}' ---")
            
            # --- Print Statistics ---
            print('\n--- Sequence Statistics -----')
            print(f'DNA sequence length: {sequence_info.dna_seq_len} nucleotides')
            print(f'GC content: {sequence_info.dna_gc_content:.2f}%')

            # --- ORFs ---
            print('\n----- Open Reading Frames (ORFs) -----')
            print(f'Number of ORFs found: {sequence_info.orfs_amount}')
            if sequence_info.orfs:
                for orf_idx, orf in enumerate(sequence_info.orfs, start=1):
                    print(f'   ORF {orf_idx}: {orf} ({len(orf)} bp, GC: {sequence_info.orfs_gc[orf_idx-1]:.2f}%)')
                
                # --- Transcription & Translation ---
                print('\n----- Transcription & Translation -----')
                print(f'Transcribed RNA from full sequence: {sequence_info.raw_rna_seq}')
                print(f'Translated Protein from full sequence: {sequence_info.protein_seq}\n')
                
                for orf_idx, (orf_rna, orf_protein) in enumerate(zip(sequence_info.orfs_rna, sequence_info.orfs_protein), start=1):
                    print(f'ORF {orf_idx} RNA: {orf_rna}')
                    print(f'ORF {orf_idx} Protein: {orf_protein}\n')
            else:
                print("No Open Reading Frames (ORFs) of minimum length found.")

            # --- Report Generation ---
            print('----- Report Generation -----')
            today = datetime.datetime.now().strftime("%Y-%m-%d")
            
            # Sanitize header for use in filename
            sanitized_header = "".join(c for c in sequence_info.header if c.isalnum() or c in ('_', '-')).rstrip()
            base_name = os.path.join(REPORTS_DIR, f"analysis_{sanitized_header}_{today}")
            
            unique_base_name = create_report.generate_report_filename(base_name)
            print(f"Report files will be based on: {unique_base_name}")

            # Create and write text report
            report_content = create_report.create_report_content(sequence_info)
            if create_report.write_report_to_file(report_content, unique_base_name):
                print(f"Text report saved: {unique_base_name}.txt")
            else:
                print("Failed to create the text report.")

            # Serialize SequenceInfo to a dictionary and write to JSON
            if create_report.write_data_to_json(sequence_info.to_dict(), unique_base_name):
                print(f"JSON report saved: {unique_base_name}.json")
            else:
                print("Failed to create the JSON report.")

    except (ValueError, FileNotFoundError) as e:
        print(f"\nAn error occurred: {e}")
    except Exception as e:
        # Catch any other unexpected errors
        print(f"\nAn unexpected error occurred: {e}")

if __name__ == "__main__":
    main()