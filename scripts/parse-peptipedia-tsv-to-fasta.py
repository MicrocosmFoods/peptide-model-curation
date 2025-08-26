import argparse
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

# Define the allowed amino acid alphabet
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

def is_valid_sequence(sequence):
    """Check if the sequence contains only valid amino acids."""
    return set(sequence).issubset(valid_amino_acids)

def parse_tsv(input_file):
    """Parse the TSV file and return a list of valid sequence records."""
    sequence_records = []
    
    with open(input_file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        next(reader)  # Skip header
        for row in reader:
            peptide_id = row[0]
            peptide_sequence = row[1]
            
            if is_valid_sequence(peptide_sequence):
                sequence_records.append((peptide_id, peptide_sequence))
            else:
                print(f"Warning: Skipping invalid sequence {peptide_sequence} for ID {peptide_id}")
    
    return sequence_records

def write_fasta(output_file, records):
    """Write sequences to a FASTA file."""
    seq_records = [
        SeqRecord(Seq(seq), id=id_, description="")
        for id_, seq in records
    ]
    with open(output_file, 'w') as fastafile:
        SeqIO.write(seq_records, fastafile, 'fasta')

def main():
    parser = argparse.ArgumentParser(description="Convert a TSV file of peptide sequences to FASTA format")
    parser.add_argument("input", help="Input TSV file with peptide ID and sequence")
    parser.add_argument("output", help="Output FASTA file")
    args = parser.parse_args()

    # Parse the TSV file and retrieve valid sequences
    sequence_records = parse_tsv(args.input)

    # Write the valid sequences to the output FASTA file
    write_fasta(args.output, sequence_records)

    print(f"FASTA file '{args.output}' created successfully with {len(sequence_records)} sequences.")

if __name__ == "__main__":
    main()