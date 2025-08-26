import argparse
import csv
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import OrderedDict

# Define the allowed amino acid alphabet
valid_amino_acids = set("ACDEFGHIKLMNPQRSTVWY")

def is_valid_sequence(sequence):
    """Check if the sequence contains only valid amino acids."""
    return set(sequence).issubset(valid_amino_acids)

def parse_tsv(input_file, metadata_file):
    """Parse the TSV files and return OrderedDicts for sequences and metadata."""
    sequence_records = OrderedDict()
    metadata_records = OrderedDict()
    
    # Parse the input TSV file
    with open(input_file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        next(reader)  # Skip header
        for row in reader:
            fmdb_id = row[0]
            peptide_sequence = row[1]
            
            if not is_valid_sequence(peptide_sequence):
                print(f"Warning: Skipping invalid sequence {peptide_sequence} for ID {fmdb_id}")
                continue
            
            if peptide_sequence not in sequence_records:
                sequence_records[peptide_sequence] = fmdb_id

    # Parse the metadata TSV file
    with open(metadata_file, 'r') as metafile:
        reader = csv.DictReader(metafile, delimiter='\t')
        for row in reader:
            fmdb_id = row['FMDB_ID']
            if fmdb_id in sequence_records.values():
                metadata_records[fmdb_id] = row

    return sequence_records, metadata_records

def write_fasta(output_file, records):
    """Write sequences to a FASTA file."""
    seq_records = [
        SeqRecord(Seq(seq), id=id_, description="")
        for seq, id_ in records.items()
    ]
    with open(output_file, 'w') as fastafile:
        SeqIO.write(seq_records, fastafile, 'fasta')

def write_metadata_tsv(output_file, metadata_records):
    """Write the filtered metadata to a new TSV file."""
    if metadata_records:
        with open(output_file, 'w', newline='') as tsvfile:
            fieldnames = metadata_records[next(iter(metadata_records))].keys()
            writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter='\t')
            writer.writeheader()
            for row in metadata_records.values():
                writer.writerow(row)

def main():
    parser = argparse.ArgumentParser(description="Convert a TSV file of peptide sequences to FASTA format and filter metadata")
    parser.add_argument("input", help="Input TSV file with FMDB_ID and Peptide_Sequence")
    parser.add_argument("metadata", help="Input metadata TSV file")
    parser.add_argument("output_fasta", help="Output FASTA file")
    parser.add_argument("output_metadata", help="Output filtered metadata TSV file")
    args = parser.parse_args()

    # Parse the TSV files and retrieve unique, valid sequences and corresponding metadata
    sequence_records, metadata_records = parse_tsv(args.input, args.metadata)

    # Write the valid sequences to the output FASTA file
    write_fasta(args.output_fasta, sequence_records)

    # Write the filtered metadata to the output TSV file
    write_metadata_tsv(args.output_metadata, metadata_records)

    print(f"FASTA file '{args.output_fasta}' created successfully with {len(sequence_records)} unique sequences.")
    print(f"Filtered metadata TSV file '{args.output_metadata}' created successfully with {len(metadata_records)} records.")

if __name__ == "__main__":
    main()
