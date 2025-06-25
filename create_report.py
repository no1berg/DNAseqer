# This module provides functions to create a report from DNA sequence analysis data.

import datetime
from email.mime import base
import json
import os
import re

# Generate a unique filename for the report
def generate_report_filename(base_name: str = 'report') -> str:
    """
    Generate a unique filename for the report by appending a counter if the file already exists.

    This function checks for the existence of both the txt and json files with the same base name, 
    to avoid overwriting existing reports. It appends a counter to the base name if necessary.
    Args:
        base_name (str): The base name for the report file (default is 'report').
    Returns:
        str: A unique filename for the report, with a timestamp and appropriate extension.
    """
    # Check if the initial name for either extension already exists
    if not os.path.exists(f"{base_name}.txt") and not os.path.exists(f"{base_name}.json"):
        return base_name
    
    #-- If the file exists, append a counter to the base name
    counter = 1
    while True:
        # Generate a new filename with the counter
        new_name = f"{base_name}_{counter}"
        # Check for both txt and json files
        if not os.path.exists(f"{new_name}.txt") and not os.path.exists(f"{new_name}.json"):
            return new_name
        counter += 1

# Create a report content from the provided data
def create_report_content(data: dict) -> str:
    """Generate a report content from the provided data.
    This function takes a dictionary of data and formats it into a string report.

    Args:
        data (dict): A dictionary containing the report data.

    Returns:
        str: The formatted report as a string.
    """
    #-- Create a header for the report
    content = [
        f"Report generated on: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        f"DNA Sequence Report\n",
        "Analysis Results:\n",
        "-" * 40 + "\n",
    ]

    #-- Data processing and handling
    if not data:
        content.append("No data available for report generation.\n")
    else:
        content.append(f"DNA Sequence: {data.get('dna_seq', 'N/A')}\n")
        content.append(f"DNA Sequence Length: {data.get('dna_seq_len', 'N/A')} nucleotides\n")
        content.append(f"GC Content: {data.get('dna_gc_content', 'N/A')}%\n")
        content.append(f"Transcribed RNA Sequence: {data.get('raw_rna_seq', 'N/A')}\n")
        content.append("-" * 40 + "\n")
        content.append(f"Number of Open Reading Frames (ORFs): {data.get('orfs_amount', 0)}\n")
        content.append("Open Reading Frames (ORFs): \n")
        for i, orf in enumerate(data.get('orfs', [])):
            content.append(f"  - ORF {i+1} ({len(orf)}bp): {orf}\n")
        #-- Add ORF cg content if available
        if 'orfs_gc' in data:
            content.append("GC Content for each ORF:\n")
            for i, gc in enumerate(data.get('orfs_gc', [])):
                content.append(f"  - ORF {i+1}: {gc:.2f}%\n")
        else:
            content.append("No GC content data for ORFs.\n")
        content.append("-" * 40 + "\n")

        #-- Add ORF RNA sequences if available
        for i, rna in enumerate(data.get('orfs_rna', [])):
            content.append(f"Transcribed RNA from ORF {i+1}: {rna}\n")
        content.append("-" * 40 + "\n")
        #-- Add ORF protein sequences if available
        for i, protein in enumerate(data.get('orfs_protein', [])):
            content.append(f"Translated Protein from ORF {i+1}: {protein}\n")
        content.append("-" * 40 + "\n")
        content.append(f"Protein Sequence: {data.get('protein_seq', 'N/A')}\n")
        # calculate the protein length
        protein_length = len(data.get('protein_seq', ''))
        content.append(f"Protein Sequence Length: {protein_length} amino acids\n")
    #-- Add a footer to the report
    content.extend([
        "-" * 40 + "\n",
        "End of Report\n"
    ])

    return ''.join(content)

def write_report_to_file(content: str, file_name) -> bool:
    """Writes the report content to a specified file.

    Args:
        file_name (str): The name of the file to write the report to (without extension).
        content (str): The content of the report to write.
    Returns:
        bool: True if the report was successfully written, False otherwise.
    """
    try:
        with open(f"{file_name}.txt", 'w') as file:
            file.write(content)
            print(f"Successfully written report to {file_name}")
        return True
    except IOError as e:
        print(f"Error writing report to file: {e}")
        return False
    
# Store the report content in a JSON file
def write_data_to_json(data: dict, file_name: str) -> bool:
    """ Writes the report data to a JSON file.
    
    Args:
        data (dict): The data to write to the JSON file.
        file_name (str): The name of the JSON file to write to.
    Returns:
        bool: True if the data was successfully written, False otherwise.
    """
    try:
        with open(f"{file_name}.json", 'w') as file:
            json.dump(data, file, indent=4)
            print(f"Successfully written data to {file_name}")
        return True
    except (IOError, TypeError) as e:
        # Catch IOError for file writing issues and TypeError for JSON serialization issues
        print(f"Error: Could not write data to JSON file: {file_name}. Error: {e}")
        return False