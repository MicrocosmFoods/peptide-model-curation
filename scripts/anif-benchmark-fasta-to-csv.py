#!/usr/bin/env python3
import argparse
from Bio import SeqIO
import csv

def convert_fasta_to_csv(input_fasta, output_csv):
    """
    Convert FASTA file to CSV with sequence and bioactivity columns.
    Bioactivity is 1 for positive sequences and 0 for negative sequences.
    """
    # Open output CSV file
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Write header
        writer.writerow(['sequence', 'bioactivity'])
        
        # Parse FASTA file and write rows
        for record in SeqIO.parse(input_fasta, "fasta"):
            sequence = str(record.seq)
            # Check if header contains "positive" or "negative"
            bioactivity = 1 if "positive" in record.description.lower() else 0
            writer.writerow([sequence, bioactivity])

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert FASTA file to CSV with sequence and bioactivity columns')
    parser.add_argument('input', help='Input FASTA file')
    parser.add_argument('output', help='Output CSV file')
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert file
    convert_fasta_to_csv(args.input, args.output)
    print(f"Converted {args.input} to {args.output}")

if __name__ == "__main__":
    main()