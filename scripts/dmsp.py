import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import numpy as np
import tensorflow as tf
from tensorflow import keras

# Argument parser to handle input and output files
parser = argparse.ArgumentParser(description='''Predicts the detectability of input peptides using a single dimension
                                                Convolutional Neural Network, based on Tensorflow 1.13.1
                                                Requirements: Tensorflow 1.13.1''')
parser.add_argument('infile', metavar='F', type=str,
                    help='File containing the peptides to be predicted, one per line (max length= 81)')
parser.add_argument('outfile', type=str,
                    help='File to save the prediction results')
args = parser.parse_args()

# Function to load and codify peptides in chunks
def load_pep_and_codify_chunk(file, max_len, start_idx, end_idx):
    aa_dict = {
        'A': 1, 'R': 2, 'N': 3, 'D': 4, 'C': 5, 'Q': 6, 'E': 7, 'G': 8, 'H': 9, 'I': 10, 'L': 11, 'K': 12,
        'M': 13, 'F': 14, 'P': 15, 'O': 16, 'S': 17, 'U': 18, 'T': 19, 'W': 20, 'Y': 21, 'V': 22
    }
    with open(file, 'r') as inf:
        lines = inf.read().splitlines()[start_idx:end_idx]
    
    pep_codes = []
    long_pep_counter = 0
    invalid_pep_counter = 0
    newLines = []
    accessions = []
    
    for line in lines:
        accession, pep = line.split('\t')
        if 'X' in pep or any(aa not in aa_dict for aa in pep):  # Check 'X' or any invalid character
            invalid_pep_counter += 1
            continue
        if len(pep) > max_len:
            long_pep_counter += 1
            continue
        current_pep = [aa_dict[aa] for aa in pep]
        pep_codes.append(current_pep)
        newLines.append(pep)
        accessions.append(accession)
    
    predict_data = keras.preprocessing.sequence.pad_sequences(pep_codes, value=0, padding='post', maxlen=max_len)
    return predict_data, long_pep_counter, invalid_pep_counter, newLines, accessions

# Load the model
try:
    print('Loading model...')
    model_2_1D = keras.models.load_model('model/model_2_1D.h5') 
except Exception as e:
    print(f'Error loading model: {e}')
    exit(1)

# Function to count the total number of lines in the input file
def count_lines_in_file(file):
    with open(file, 'r') as f:
        return sum(1 for _ in f)

# Main processing loop for chunking
chunk_size = 100000  # Adjust this based on your available memory
total_peptides = count_lines_in_file(args.infile)
start_idx = 0

total_long_pep_counter = 0
total_invalid_pep_counter = 0
total_processed = 0

# Open output file and write the header
with open(args.outfile, 'w') as outf:
    outf.write('Header\tPeptide\tProb\tDetectability\n')
    
    while start_idx < total_peptides:
        end_idx = min(start_idx + chunk_size, total_peptides)
        
        # Load and process a chunk of peptides
        predict_data, long_pep_counter, invalid_pep_counter, lines, accessions = load_pep_and_codify_chunk(args.infile, 81, start_idx, end_idx)
        
        if len(lines) == 0:  # Skip empty chunks (in case all peptides were invalid or too long)
            start_idx += chunk_size
            continue
        
        # Make predictions on the chunk
        model_2_1D_pred = model_2_1D.predict(predict_data)
        
        # Process and save predictions
        model_2_1D_pred = np.hstack((np.array(lines).reshape(len(lines), 1), model_2_1D_pred)).tolist()
        
        Pred_output = []
        for i, pred in enumerate(model_2_1D_pred):
            if float(pred[1]) > 0.5:
                Pred_output.append([accessions[i], pred[0], str(1 - float(pred[1])), '0'])
            else:
                Pred_output.append([accessions[i], pred[0], str(1 - float(pred[1])), '1'])
        
        # Append predictions to the output file
        with open(args.outfile, 'a') as outf:
            outf.writelines('\t'.join(i) + '\n' for i in Pred_output)
        
        # Update counters
        total_long_pep_counter += long_pep_counter
        total_invalid_pep_counter += invalid_pep_counter
        total_processed += len(lines)

        # Calculate percentages
        total_attempted = total_processed + total_long_pep_counter + total_invalid_pep_counter
        loaded_percentage = (len(lines) / total_attempted) * 100
        long_pep_percentage = (long_pep_counter / total_attempted) * 100
        invalid_pep_percentage = (invalid_pep_counter / total_attempted) * 100

        print(f'Processed chunk {start_idx} to {end_idx}:')
        print(f'  - Peptides loaded: {len(lines)} ({loaded_percentage:.2f}%)')
        print(f'  - Long peptides skipped: {long_pep_counter} ({long_pep_percentage:.2f}%)')
        print(f'  - Invalid peptides skipped: {invalid_pep_counter} ({invalid_pep_percentage:.2f}%)')
        
        start_idx += chunk_size

print('Finished processing all peptides.')
print(f'Total peptides processed: {total_processed}')
print(f'Total long peptides skipped: {total_long_pep_counter}')
print(f'Total invalid peptides skipped: {total_invalid_pep_counter}')

