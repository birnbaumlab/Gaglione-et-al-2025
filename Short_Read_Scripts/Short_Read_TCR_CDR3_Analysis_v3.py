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
import re


def generate_reverse_reads(reverse_file):
    """
    Generator function to yield sequences from a FASTQ file with reverse reads (Read 2).
    In our sequencing format, the CDR3 is entirely contained in the reverse read.
    This function takes each reverse read, takes reverse complement, and yields the sequence in a generator
    """
    if reverse_file.endswith(".fastq"):
        #print(f"Reading file: {file_path}")
        with open(reverse_file, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                yield str(record.seq.reverse_complement())
    else:
        raise ValueError('Invalid file type')


def generate_batches(generator, batch_size):
    """
    Generator helper function to yield batches of sequences from a generator.
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
    Helper function to create a printed quality report for a batch of data.
    """
    total = len(data)
    print('-- BATCH REPORT --')
    good = sum([(d != False and d != None) for d in data])
    print('> {} ({}%) [{}]'.format(good,round(100*good/total,1),'Matched'))
    print('-------')
    print('{} entries in batch'.format(total))


def create_barcode_dict(params):
    """
    Create a dictionary to map CDR3 sequences to their corresponding index/identity for counting.
    The dictionary is created from a reference file containing CDR3 sequences and their corresponding indices.
    The function also generates SNP variants of the CDR3 sequences and adds them to the dictionary.
    The SNP variants are created by replacing each nucleotide in the CDR3 sequence with all possible nucleotides (A, T, G, C)
    and labeled as either '_silent' or '_mutant' based on whether the amino acid sequence remains the same or not.
    """

    cdr3_barcode_dict = {}
    ref_file = pd.read_csv(params['Reference_File'])
    for index, row in ref_file.iterrows():
        seq = row.Sequence.upper()
        #use regex match to find index of constant flank
        flank_pos = seq.find(params['Constant_Flank'])
        if flank_pos == -1:
            #print(seq)
            #raise ValueError('Error in making barcode dictionary.')
            continue
        cdr3_barcode = seq[flank_pos - params['Barcode_Offset'] - params['Window_Size']: flank_pos - params['Barcode_Offset']]
        cdr3_barcode_dict[cdr3_barcode] = str(index + 1)

        #create SNP variants of the barcode
        assert params['Window_Size'] % 3 == 0
        cdr3_barcode_aminos = Seq(cdr3_barcode).translate()
        for i in range(len(cdr3_barcode)):
            for nt in ['A', 'T', 'G', 'C']:
                if cdr3_barcode[i] == nt:
                    continue
                cdr3_barcode_snp = cdr3_barcode[:i] + nt + cdr3_barcode[i+1:]
                if Seq(cdr3_barcode_snp).translate() == cdr3_barcode_aminos:
                    snp_label = '_silent'
                else:
                    snp_label = '_mutant'

                if cdr3_barcode_snp not in cdr3_barcode_dict:
                    cdr3_barcode_dict[cdr3_barcode_snp] = str(index + 1) + snp_label
    
    return cdr3_barcode_dict


def process_read_batch(sequences, params, cdr3_barcode_dict, verbose=True):
    """
    Process a batch of sequences to extract CDR3 sequences and match them with the barcode dictionary.
    """

    start_time = time.perf_counter()
    processed_batch = []

    for seq in sequences:
        seq = seq[:-1*(params['UMI_Size'])]
        flank_pos = len(seq)-1*len(params['Constant_Flank'])
        result = None
        if seq[flank_pos:] != params['Constant_Flank']:
            #print(seq[flank_pos:], params['Constant_Flank'])
            #print(seq)
            continue
        else:
            try:
                barcode = seq[flank_pos - params['Barcode_Offset'] - params['Window_Size']: flank_pos - params['Barcode_Offset']]
                result = cdr3_barcode_dict[barcode]
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
    parser.add_argument('params', type=str, help='Path to params file')
    parser.add_argument('input_fastq_read2', type=str, help='Path to FASTQ file containing read 2 (reverse reads)')
    parser.add_argument('output_handle', type=str, help='Handle name and directory to output file')
    args = parser.parse_args()

    with open(args.params, 'r') as f:
        params = json.load(f)

    # Initialize generators
    #all_reads = generate_forward_reads(args.input_fastq_read1)
    all_reads = generate_reverse_reads(args.input_fastq_read2)
    read_batches = generate_batches(all_reads, params['Batch_Size'])
    cdr3_barcode_dict = create_barcode_dict(params)
    cdr3_counts_dict = {}

    with open(f'{args.output_handle}.txt', "a") as outfile:
        outfile.write('CDR3_Match\n')
        
        # Use joblib with return_as="generator"
        results_gen = Parallel(params['N_Processes'], return_as="generator")(delayed(process_read_batch)(batch, params, cdr3_barcode_dict, verbose=True) for batch in read_batches)
        #save 

        # Write each processed chunk to the file
        for chunk in tqdm(results_gen):
            outfile.writelines(str(item) + '\n' for item in chunk)
            cdr3_counts_dict = append_chunk_to_counts_dict(cdr3_counts_dict, chunk)

            with open(f'{args.output_handle}.json', 'w') as json_outfile:
                json.dump(cdr3_counts_dict, json_outfile)