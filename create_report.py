## This module provides functions to create a report from DNA sequence analysis data.

import datetime
import json
import os
from typing import Dict, Any

# Import SequenceInfo for type hinting to avoid circular dependencies
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from main import SequenceInfo


def generate_report_filename(base_name: str = 'report') -> str:
    """
    Generates a unique filename by appending a counter if a .txt or .json file exists.
    """
    if not os.path.exists(f"{base_name}.txt") and not os.path.exists(f"{base_name}.json"):
        return base_name
    
    counter = 1
    while True:
        new_name = f"{base_name}_{counter}"
        if not os.path.exists(f"{new_name}.txt") and not os.path.exists(f"{new_name}.json"):
            return new_name
        counter += 1


def create_report_content(data: 'SequenceInfo') -> str:
    """
    Generates a formatted string report from a SequenceInfo object.

    Args:
        data (SequenceInfo): An object containing the analysis data.

    Returns:
        str: The formatted report as a string.
    """
    # Access attributes directly from the object, not with .get()
    content = [
        f"DNA Sequence Report\n",
        f"Report generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        "-" * 50 + "\n",
        "--- SEQUENCE STATISTICS ---\n",
        f"Original DNA Sequence: {data.dna_seq}\n",
        f"Sequence Length: {data.dna_seq_len} nucleotides\n",
        f"GC Content: {data.dna_gc_content:.2f}%\n",
        f"Full Transcribed RNA: {data.raw_rna_seq}\n",
        f"Translated Protein (from first start codon): {data.protein_seq}\n",
        f"Protein Length: {len(data.protein_seq)} amino acids\n",
        "-" * 50 + "\n",
        "--- OPEN READING FRAMES (ORFs) ---\n",
        f"Number of ORFs found (min. 50bp): {data.orfs_amount}\n\n"
    ]

    if not data.orfs:
        content.append("No Open Reading Frames (ORFs) found.\n")
    else:
        for i, orf in enumerate(data.orfs, 1):
            content.append(f"ORF {i}:\n")
            content.append(f"  - DNA ({len(orf)} bp): {orf}\n")
            # Correctly access the corresponding GC content
            if i <= len(data.orfs_gc):
                content.append(f"  - GC Content: {data.orfs_gc[i-1]:.2f}%\n")
            # Correctly access the corresponding RNA and Protein
            if i <= len(data.orfs_rna):
                content.append(f"  - RNA: {data.orfs_rna[i-1]}\n")
            if i <= len(data.orfs_protein):
                content.append(f"  - Protein: {data.orfs_protein[i-1]}\n\n")

    content.extend([
        "-" * 50 + "\n",
        "End of Report\n"
    ])

    return ''.join(content)


def write_report_to_file(content: str, file_name: str) -> bool:
    """
    Writes the report content to a specified text file.
    """
    try:
        with open(f"{file_name}.txt", 'w') as file:
            file.write(content)
        print(f"Successfully written report to {file_name}.txt")
        return True
    except IOError as e:
        print(f"Error writing report to file: {e}")
        return False


def write_data_to_json(data: Dict[str, Any], file_name: str) -> bool:
    """
    Writes dictionary data to a JSON file.
    """
    try:
        with open(f"{file_name}.json", 'w') as file:
            json.dump(data, file, indent=4)
        print(f"Successfully written data to {file_name}.json")
        return True
    except (IOError, TypeError) as e:
        print(f"Error: Could not write data to JSON file '{file_name}.json'. Error: {e}")
        return False