#!/usr/bin/env python
# -*- coding: utf-8 -*-

import gzip
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq, reverse_complement
import os
from tqdm import tqdm
import numpy as np
from joblib import Parallel, delayed
import argparse
import time
import json


def generate_sequences(dir_path):
    """
    Iterate through all fastq.gz files in directory 'dir_path' and yield sequences through a generator.
    """
    for file_path in os.listdir(dir_path):
        if file_path.endswith(".fastq.gz"):
            with gzip.open(os.path.join(dir_path, file_path), "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    yield str(record.seq)


def generate_batches(generator, batch_size):
    """
    Yield batches of size 'batch_size' from the generator.
    """

    batch = []
    for item in generator:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def create_cdr3a_cdr3b_barcode_dicts(params):
    """
    Create a dictionary to map CDR3 sequences to their corresponding index/identity for counting.
    The dictionary is created from a reference file containing CDR3 sequences and their corresponding indices.
    The function also generates SNP variants of the CDR3 sequences and adds them to the dictionary.
    The SNP variants are created by replacing each nucleotide in the CDR3 sequence with all possible nucleotides (A, T, G, C)
    and labeled as either '_silent' or '_mutant' based on whether the amino acid sequence remains the same or not.
    """
    def create_barcode_dict(bc_type):
        cdr3_barcode_dict = {}
        ref_file = pd.read_csv(params['Reference_File'])
        for index, row in ref_file.iterrows():
            seq = row.Sequence

            #use regex match to find index of constant flank
            flank_pos = np.max([seq.find(i) for i in params[f'{bc_type}_Barcode_Constant_Flanks']])
            if flank_pos == -1:
                #print(seq)
                #raise ValueError('Error in making barcode dictionary.')
                continue
            cdr3_barcode = seq[flank_pos + params[f'{bc_type}_Barcode_Offset'] : flank_pos + params[f'{bc_type}_Barcode_Offset'] + params[f'{bc_type}_Barcode_Size']]
            cdr3_barcode = cdr3_barcode.upper()
            cdr3_barcode_dict[cdr3_barcode] = index + 1

        n_unique = np.unique(np.array(list(cdr3_barcode_dict.keys())))
        if n_unique.size != len(cdr3_barcode_dict):
            print(f'Warning: Duplicate {bc_type} barcodes found. This may cause issues with downstream analysis.')
            print(f'Number of barcode clashes: {len(cdr3_barcode_dict) - n_unique.size}')
            print(f'Recommend increasing the parameter "{bc_type}_Barcode_Size" in the params file.')

        print(f'Created {bc_type} barcode dictionary.')
        return cdr3_barcode_dict

    return create_barcode_dict('CDR3A'), create_barcode_dict('CDR3B')


def quality_report(data):
    """
    Print a quality report of the processed data.
    """

    total = len(data)
    print('-- BATCH REPORT --')
    for key in data[0].keys():
        good = sum([(d[key] != False and d[key] != None) for d in data])
        print('> {} ({}%) [{}]'.format(good,round(100*good/total,1),key))
    print('-------')
    print('{} entries in batch'.format(total))


def format_read_row(read_row_dict):
    """
    Format a single annotated read dictionary into a CSV row.
    """

    return str(read_row_dict['Exceeds_Size_Cutoff']) + ',' + \
           str(read_row_dict['Contains_TRAC']) + ',' + \
           str(read_row_dict['Contains_TRBC']) + ',' + \
           str(read_row_dict['Matched_TRAV']) + ',' + \
           str(read_row_dict['Matched_TRBV']) + ',' + \
           str(read_row_dict['Matched_CDR3A_BC']) + ',' + \
           str(read_row_dict['Matched_CDR3B_BC']) + ',' + \
           str(read_row_dict['Flagged']) + ',' + \
           str(read_row_dict['Sequence']) + '\n'


def process_read_batch(batch, params, cdr3_barcodes, verbose = True):
    """
    Run full annotation pipeline on a batch of reads, and return a list of dictionaries with the annotated results.
    """

    processed_batch = []

    # add a timer
    start_time = time.perf_counter()

    # Read all the dictionaries 1 time, not inside of each read iteration
    all_TCR_list = pd.read_csv(params['Reference_File'])
    trbv_dict_unique = pd.read_csv(params['TRBV_dict_unique'])
    second_trbv_dict = pd.read_csv(params['TRBV_dict'])    
    trav_dict_unique = pd.read_csv(params['TRAV_dict_unique'])
    second_trav_dict = pd.read_csv(params['TRAV_dict'])

    cdr3a_dict, cdr3b_dict = cdr3_barcodes
    cdr3a_len, cdr3b_len = (params['CDR3A_Barcode_Size'], params['CDR3B_Barcode_Size'])
    
    for i,read in enumerate(batch):

        result = {'Exceeds_Size_Cutoff': len(read) >= params['Size_Cutoff'],
                  'Contains_TRAC': False,
                  'Contains_TRBC': False,
                  'Extracted_CDR3A': None,
                  'Extracted_CDR3B': None,
                  'Matched_TRBV': None,
                  'Matched_TRAV': None,
                  'Matched_CDR3B_BC': None,
                  'Matched_CDR3A_BC': None,
                  'Flagged': False,
                  'Sequence': None}
        
        #1. Filter by size
        if result['Exceeds_Size_Cutoff']:

            reverse_read = reverse_complement(read)

            #2. Check for TRAC substrings first (outermost component)
            for substring, index in params['TRAC_substrings']:
                trac_pos = read.find(substring)
                if trac_pos == -1:
                    trac_pos = reverse_read.find(substring)
                    if trac_pos != -1:
                        read = reverse_read

                if trac_pos != -1:
                    result['Contains_TRAC'] = True
                    result['Extracted_CDR3A'] = read[trac_pos - index - params['CDR3_search_window']: trac_pos - index]
                    break

            
            #3. Check for TRBC substrings - reverse reads are already checked above
            for substring, index in params['TRBC_substrings']:
                trbc_pos = read.find(substring)
                if trbc_pos != -1:
                    result['Contains_TRBC'] = True
                    result['Extracted_CDR3B'] = read[trbc_pos - index - params['CDR3_search_window']: trbc_pos - index]
                    break


            if result['Contains_TRAC'] and result['Contains_TRBC']:

                read = read[:trac_pos + 30] #narrows down the search space for TRAV and TRBV

                #4. Check for TRBV substrings
                for _, row_unique in trbv_dict_unique.iterrows():
                    trbv_allele_unique = row_unique['TRBV']
                    barcode_unique = row_unique['TRBV_barcode']
                    barcode_pos = read.find(barcode_unique)
                    if barcode_pos != -1:
                        second_TRBV_barcode = second_trbv_dict[second_trbv_dict.TRBV == trbv_allele_unique]['TRBV_barcode'].values[0] #2nd TRBV barcode to confirm the TRBV label 
                        if second_TRBV_barcode in read:
                            result['Matched_TRBV'] = trbv_allele_unique
                            break    

                #5. Check for TRAV substrings
                for _, row_unique in trav_dict_unique.iterrows():
                    trav_allele_unique = row_unique['TRAV']
                    barcode_unique = row_unique['TRAV_barcode']
                    barcode_pos = read.find(barcode_unique)
                    if barcode_pos != -1:
                        second_TRAV_barcode = second_trav_dict[second_trav_dict.TRAV == trav_allele_unique]['TRAV_barcode'].values[0] #2nd TRAV barcode to confirm the TRAV label 
                        if second_TRAV_barcode in read:
                            result['Matched_TRAV'] = trav_allele_unique
                            break

                
                #6a. Check for CDR3A barcodes:
                for i in range(0,len(result['Extracted_CDR3A'])):
                    cdr3a_fragment = result['Extracted_CDR3A'][i:i+cdr3a_len]
                    try:
                        result['Matched_CDR3A_BC'] = cdr3a_dict[cdr3a_fragment]
                        break
                    except KeyError: continue
                        
                #6b. Check for CDR3B barcodes:
                for i in range(0,len(result['Extracted_CDR3B'])):
                    cdr3b_fragment = result['Extracted_CDR3B'][i:i+cdr3b_len]
                    try:
                        result['Matched_CDR3B_BC'] = cdr3b_dict[cdr3b_fragment]
                        break
                    except KeyError: continue     
                           

                #7. Store the original sequence if CDR3 TCR numbers don't match or if there is a mismatched TCR
                if result['Matched_CDR3A_BC'] is not None and result['Matched_CDR3B_BC'] is not None:
                    if result['Matched_CDR3A_BC'] != result['Matched_CDR3B_BC']:
                        result['Sequence'] = read
                        result['Flagged'] = True

                    elif result['Matched_TRBV'] is not None and result['Matched_TRAV'] is not None:
                        proper_TRBV_pairing = all_TCR_list.loc[result['Matched_CDR3B_BC']-1, 'TRBV'] == result['Matched_TRBV']
                        proper_TRAV_pairing = all_TCR_list.loc[result['Matched_CDR3B_BC']-1, 'TRAV'] == result['Matched_TRAV']

                        if not proper_TRBV_pairing or not proper_TRAV_pairing:
                            result['Sequence'] = read
                            result['Flagged'] = True
                
        processed_batch.append(result)
    
    # run quality report
    if verbose:
        quality_report(processed_batch)

        execution_time = time.perf_counter() - start_time
        print('Batch completed in: {} s'.format(round(execution_time,2)))
    
    return processed_batch



if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('params', type=str, help='Path to params file')
    parser.add_argument('input_dir', type=str, help='Path to input dir')
    parser.add_argument('output_csv', type=str, help='Path to output CSV file')
    args = parser.parse_args()

    with open(args.params, 'r') as f:
        params = json.load(f)

    # Initialize generators
    all_reads = generate_sequences(args.input_dir)
    read_batches = generate_batches(all_reads, params['batch_size'])

    # initialize barcodes
    cdr3_barcode_output = create_cdr3a_cdr3b_barcode_dicts(params)

    with open(args.output_csv, "a") as outfile:
        outfile.write('Exceeds_Size_Cutoff,Contains_TRAC,Contains_TRBC,Matched_TRAV,Matched_TRBV,Matched_CDR3A_BC,Matched_CDR3B_BC,Flagged,Sequence\n')
        
        # Use joblib with return_as="generator"
        results_gen = Parallel(n_jobs=params['n_processes'], return_as="generator")(delayed(process_read_batch)(batch, params, cdr3_barcode_output) for batch in read_batches)
        #save 

        # Write each processed chunk to the file
        for chunk in tqdm(results_gen):
            outfile.writelines(format_read_row(item) for item in chunk)