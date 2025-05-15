#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
from tqdm import tqdm
import argparse


def generate_sequences(file_path):
    """
    Iterate through a fastq file and yield sequences one by one as a generator.
    param file_path: Path to the fastq file
    return: A generator that yields sequences from the fastq file
    """
    if file_path.endswith(".fastq"):
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                yield str(record.seq)


def annotate_TRAV_or_TRBV(generator, v_dict_unique_path, trav_or_trbv, origin_flank):
    """
    Annotate TRAV or TRBV sequences from the generator.
    params:
        generator: A generator that yields sequences
        v_dict_unique_path: Path to the unique V gene dictionary
        trav_or_trbv: Specify TRAV or TRBV
        origin_flank: The origin flank sequence to orient the reads
        search_cutoff: The cutoff length for searching for the barcode
    return: 
        A generator that yields matched V alleles for each read, or None if no match is found
    """
    v_dict_unique = pd.read_csv(v_dict_unique_path)  # unique substring list

    for sequence in tqdm(generator):
        
        # If the reverse complement is found, take the reverse complement of the sequence
        if str(Seq(origin_flank).reverse_complement()) in sequence:
            sequence = str(Seq(sequence).reverse_complement())
        
        matched_v = None
        #TRAV/TRBV ANALYSIS
        # Check each TRAV/TRBV barcode
        for _, row_unique in v_dict_unique.iterrows():
            v_allele_unique = row_unique[trav_or_trbv]
            barcode_unique = row_unique[f'{trav_or_trbv}_barcode']
            
            if barcode_unique in sequence:
                matched_v = v_allele_unique
                break
        
        #if matching failed, move to the next sequence.
        if matched_v is None:
            continue

        yield matched_v


def make_vector_plot(read_list, title):
    """
    Create a vector plot from the provided read list.
    params:
        read_list: List of matched V alleles
        title: Title for the plot and output files
    return: 
        None, but saves the plot and CSV file used to make the plot
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    vector_count = pd.Series(read_list).value_counts().sort_values()
    plt.bar(vector_count.index, vector_count)
    plt.xticks(rotation=90, fontsize=10)
    plt.xlabel(title)
    plt.ylabel('Read Count')
    plt.title(f'{title} vector counts')
    plt.savefig(f'{title}_vector_counts.png', bbox_inches='tight')
    vector_count.to_csv(f'{title}_vector_counts.csv', header=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Analyze FASTQ data for the initial vector pools to confirm vector pool quality.')
    parser.add_argument('--fastq', type=str, required=True, help='Path to the fastq file')
    parser.add_argument('--v_dict_unique', type=str, required=True, help='Path to the unique V gene dictionary. Use TRBV barcodes for Step 1 and TRAV barcodes for Step 2.')
    parser.add_argument('--trav_or_trbv', type=str, required=True, choices=['TRAV', 'TRBV'], help='Specify TRAV or TRBV. Step 1 Vectors are TRBV, and Step 2 Vectors are TRAV.')
    parser.add_argument('--title', type=str, required=True, help='Title to place at top of plot and to save as')
    args = parser.parse_args()

    # Generate sequences from the fastq file
    generator = generate_sequences(args.fastq)

    origin_flank_options = {'TRAV':'TGCTGAAGCAGGCCGGTG', 'TRBV':'ATTGTTAGGTAATCGTCA'}
    origin_flank = origin_flank_options[args.trav_or_trbv]
    # Annotate TRAV or TRBV
    read_list = list(annotate_TRAV_or_TRBV(generator, args.v_dict_unique, args.trav_or_trbv, origin_flank))
    
    # Create vector plot
    make_vector_plot(read_list, args.title)
    