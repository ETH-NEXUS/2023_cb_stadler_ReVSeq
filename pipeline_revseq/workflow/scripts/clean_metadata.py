######
## Script name: clean_metadata.py
## Date: May 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the metadata and
##        removes the empty lines for controls and mouthwash samples
######

import argparse
import pandas as pd

# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='Fetch and merge multi-lane FASTA', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, help='The metadata file to clean')
    parser.add_argument('--pos', required=False, default="KoPos", type=str, help='The string to detect positive controls')
    parser.add_argument('--neg', required=False, default="KoNeg", type=str, help='The string to detect negative controls')
    parser.add_argument('--mouthwash', required=False, default="MouthWash", type=str, help='The string to detect mouthwash samples')
    parser.add_argument('--output', required=True, type=str, help='The output directory and filename')


    args = parser.parse_args()

    metadata = pd.read_table(args.input, dtype=str, sep=";")
    metadata = metadata[~metadata.ethid.str.contains(args.pos)]
    metadata = metadata[~metadata.ethid.str.contains(args.neg)]
    metadata = metadata[~metadata.ethid.str.contains(args.mouthwash)]
    
    metadata.to_csv(args.output, sep=";", index=False)

