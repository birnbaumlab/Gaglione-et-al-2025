#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
from tqdm import tqdm
import argparse
import time
import json
from joblib import Parallel, delayed 


def generate_reads(merged_file):
    """
    Generator function to yield sequences from a FASTQ file with merged forward and reverse reads.
    """

    if merged_file.endswith(".fastq"):
        #print(f"Reading file: {file_path}")
        with open(merged_file, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                yield str(record.seq)
    else:
        raise ValueError('Invalid file type')
    

def generate_batches(generator, batch_size):
    """
    Generator function to yield batches of sequences from a generator.
    """

    batch = []
    for item in generator:
        batch.append(item)
        if len(batch) == batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def quality_report(data):
    """
    Generate a quality report for a batch of data.
    """
    
    total = len(data)
    print('-- BATCH REPORT --')
    good = sum([(d != False and d != None) for d in data])
    print('> {} ({}%) [{}]'.format(good,round(100*good/total,1),'Matched'))
    print('-------')
    print('{} entries in batch'.format(total))


def create_barcode_dicts(ref_file_path, trbv_trac_frags, trbc_trav_frags, mode):
    """
    Create a dictionary to map CDR3 sequences to their corresponding index/identity for counting.
    """

    cdr3_barcode_dict = {}
    ref_file = pd.read_csv(ref_file_path)
    for index, row in tqdm(ref_file.iterrows()):
        if row.Pool != mode:
            continue

        cdr3_barcode = row.Sequence.upper()
        cdr3_barcode_dict[cdr3_barcode] = str(index + 1)

        trbv_flank = trbv_trac_frags[trbv_trac_frags.trbv == row.TRBV]['trbv_right_frag'].values[0]
        trbc_flank = trbc_trav_frags[trbc_trav_frags.trav == row.TRAV]['trbc_left_frag'].values[0]
        trav_flank = trbc_trav_frags[trbc_trav_frags.trav == row.TRAV]['trav_right_frag'].values[0]
        trac_flank = trbv_trac_frags[trbv_trac_frags.trbv == row.TRBV]['trac_left_frag'].values[0]
        
        cdr3b_start = row.Sequence.find(trbv_flank) + len(trbv_flank)
        cdr3b_end = row.Sequence.find(trbc_flank)

        cdr3a_start = row.Sequence.rfind(trav_flank) + len(trav_flank)
        cdr3a_end = row.Sequence.rfind(trac_flank)

        cdr3b_aminos = Seq(cdr3_barcode[cdr3b_start : cdr3b_end]).translate()
        cdr3a_aminos = Seq(cdr3_barcode[cdr3a_start : cdr3a_end ]).translate()
        
        try:
            assert (cdr3b_end - cdr3b_start) % 3 == 0 and (cdr3a_end - cdr3a_start) % 3 == 0
        except AssertionError:
            print(cdr3b_aminos, cdr3a_aminos)
            print(cdr3b_end - cdr3b_start, cdr3a_end - cdr3a_start)
            print(cdr3b_start, cdr3b_end, cdr3a_start, cdr3a_end)
            print(row)
            print(row.Sequence)
            print(trav_flank, trac_flank)
            
        for i in range(len(cdr3_barcode)):
            for nt in ['A', 'T', 'G', 'C']:
                if cdr3_barcode[i] == nt:
                    continue

                cdr3_barcode_snp = cdr3_barcode[:i] + nt + cdr3_barcode[i+1:]
                
                if Seq(cdr3_barcode_snp[cdr3b_start : cdr3b_end]).translate() == cdr3b_aminos and \
                   Seq(cdr3_barcode_snp[cdr3a_start : cdr3a_end]).translate() == cdr3a_aminos:
                    snp_label = '_silent'
                else:
                    snp_label = '_mutant'

                if cdr3_barcode_snp not in cdr3_barcode_dict:
                    cdr3_barcode_dict[cdr3_barcode_snp] = str(index + 1) + snp_label

    return cdr3_barcode_dict


def process_read_batch(sequences, cdr3_barcode_dict, verbose=True):
    """
    Process a batch of sequences to find matches in the CDR3 barcode dictionary.
    """

    start_time = time.perf_counter()
    processed_batch = []

    for seq in sequences:

        result = None
        try:
            result = cdr3_barcode_dict[seq[10:-8]]
        except KeyError:
            pass
        processed_batch.append(result)

    # quality report
    if verbose:
        quality_report(processed_batch)
        execution_time = time.perf_counter() - start_time
        print('Batch completed in: {} s'.format(round(execution_time,2)))
    
    return processed_batch


def append_chunk_to_counts_dict(counts_dict, chunk):
    """
    Convert a chunk of results into a dictionary for counting occurrences.
    This function takes a chunk of results and updates the counts_dict with the occurrences of each unique result.
    """

    for result in chunk:
        if result is None:
            continue
        
        result = result.split('_')[0]
        if result in counts_dict:
            counts_dict[result] += 1
        else:
            counts_dict[result] = 1
    return counts_dict


if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('input_fastq_merged', type=str, help='Path to FASTQ file containing merged reads')
    parser.add_argument('reference_CSV', type=str, help='Path to reference CSV file')
    parser.add_argument('output_handle', type=str, help='Handle name and directory to output file')
    parser.add_argument('mode', type=str, help='Pool A or Pool B - enter the letter "A" or "B".')
    args = parser.parse_args()

    # Initialize generators
    #all_reads = generate_forward_reads(args.input_fastq_read1)
    all_reads = generate_reads(args.input_fastq_merged)
    read_batches = generate_batches(all_reads, 200000)

    trbv_trac_frags = pd.read_csv('~/data/Short_Read/references/TRBV_TRAC_OH_frags.csv')
    trbc_trav_frags = pd.read_csv('~/data/Short_Read/references/TRBC_TRAV_OH_frags.csv')
    print('Starting Barcode dict creation...')
    cdr3_barcode_dict = create_barcode_dicts(args.reference_CSV, trbv_trac_frags, trbc_trav_frags, args.mode)
    print('Created Barcode Dict')
    print(len(cdr3_barcode_dict))
    cdr3_counts_dict = {}

    with open(f'{args.output_handle}.txt', "a") as outfile:
        outfile.write('CDR3_Match\n')
        
        # Use joblib with return_as="generator"
        results_gen = Parallel(n_jobs=1, return_as="generator")(delayed(process_read_batch)(batch, cdr3_barcode_dict, verbose=True) for batch in read_batches)
        #save 

        # Write each processed chunk to the file
        for chunk in tqdm(results_gen):
            outfile.writelines(str(item) + '\n' for item in chunk)
            cdr3_counts_dict = append_chunk_to_counts_dict(cdr3_counts_dict, chunk)

            with open(f'{args.output_handle}.json', 'w') as json_outfile:
                json.dump(cdr3_counts_dict, json_outfile)