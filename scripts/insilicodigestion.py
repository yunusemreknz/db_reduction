import argparse
from pyteomics import parser
from Bio import SeqIO
import os

def main():
    arg_parser = argparse.ArgumentParser(description="Perform in silico digestion on a FASTA file.")
    arg_parser.add_argument('fasta_file', type=str, help='Path to the input FASTA file')
    arg_parser.add_argument('output_file', type=str, help='Path to the output file')
    args = arg_parser.parse_args()

    fasta_file = args.fasta_file
    output_file_name = args.output_file

    # Digestion startet hier --> ungünstiges Header speichern --> CHANGE!?
    with open(output_file_name, 'w') as output_file:
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                header = record.description
                sequence = str(record.seq)
                indexed_peptides = parser.xcleave(sequence, parser.expasy_rules['trypsin'], missed_cleavages=2) # Standard trypsin digestion change if needed
                for _, peptide in indexed_peptides:
                    output_file.write(f"{header}\t{peptide}\n")

    print(f"Digested peptides with headers have been saved to '{output_file_name}'")

if __name__ == '__main__':
    main()
