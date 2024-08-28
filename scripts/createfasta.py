import argparse
import os

def wrap_sequence(sequence, line_length=60):
    """Wrap the sequence into lines of specified length."""
    return '\n'.join(sequence[i:i + line_length] for i in range(0, len(sequence), line_length))

def create_fasta_from_detected_peptides(predictions_file, output_fasta_file):
    peptide_dict = {}

    with open(predictions_file, 'r') as infile:
        infile.readline()  #skip header 

        for line in infile:
            header, peptide, prob, detectability = line.strip().split('\t')

            if detectability == '1': # filter for label=1
                if header not in peptide_dict:
                    peptide_dict[header] = peptide
                else:
                    peptide_dict[header] += peptide

    with open(output_fasta_file, 'w') as outfile:
        for i, (header, peptides) in enumerate(peptide_dict.items(), start=1):
            peptide_name = header
            
            # Save Peptides with header
            outfile.write(f">{peptide_name}\n")
            wrapped_sequence = wrap_sequence(peptides, 60)
            outfile.write(f"{wrapped_sequence}\n")

def main():
    parser = argparse.ArgumentParser(description="Generate a FASTA file from detected peptides.")
    parser.add_argument('predictions_file', type=str, help='Path to the input predictions file')
    parser.add_argument('output_fasta_file', type=str, help='Path to the output FASTA file')
    
    args = parser.parse_args()

    create_fasta_from_detected_peptides(args.predictions_file, args.output_fasta_file)
    print(f"Filtered FASTA file has been created: {args.output_fasta_file}")

if __name__ == '__main__':
    main()
