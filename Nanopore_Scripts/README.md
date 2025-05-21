# Nanopore Analysis code guide

The nanopore analysis script can be run with the following command:
`./Nanopore_analysis_v6.py /path/to/params.json /path/to/nanopore/fastq/dir path/to/put/output.csv`

All parameters for nanopore TCR annotation are found in the `params.json` file. Example params files are provided in the references folder. Here is a summary of all the parameters that can be modified in this params file:

- `n_processes`: Number of parallel processes to initiate for read annotation
- `batch_size`: Number of reads to process per parallel batch at a time
- `Reference_File`: MOST IMPORTANT! Update this to be the CDR3 oligo CSV spreadsheet created by the original TCRAFT_generate code
- `CDR3[A/B]_Barcode_Constant_Flanks`: Search flanks to find the CDR3 portions in the CDR3 oligos
- `CDR3[A/B]_Barcode_Offset`: Number of nucleotides to offset from search flank to create the unique CDR3 barcodes
- `CDR3[A/B]_Barcode_Size`: Length of the CDR3 barcodes. NOTE: You might need to increase this if the script gives a warning saying that not all barcodes are unique.
- `Size_Cutoff`: Any read smaller than this length is discarded from analysis.

We recommend keeping these parameters constant:
- `CDR3_Search_Window`: Window of the read to subset when trying to assign a CDR3 barcode to that read.
- `TRAV_dict`: Dictionary of TRAV barcodes for TRAV annotation
- `TRBV_dict`: Dictionary of TRBV barcodes for TRBV annotation
- `TRAV_dict_unique`: Longer, higher fidelity TRAV barcode set for cross-checking annotation.
- `TRBV_dict_unique`: Longer, higher fidelity TRBV barcode set for cross-checking annotation.
- `TRAC_substrings`: Search flanks and positional offsets to find the start of the TRAC
- `TRBC_substrings`: Search flanks and positional offsets to find the start of the TRBC