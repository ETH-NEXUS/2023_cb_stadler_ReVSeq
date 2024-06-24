######
## Script name: cds_consensus.py
## Date: June 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads a consensus
##        and a bed file listing all the CDS regions and outputs
##        a FASTA consensus file including only the listed regions
######

import pandas as pd
import argparse
from Bio import SeqIO


## Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='Count the number of Ns in a given consensus sequence', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, help='Input consensus sequence')
    parser.add_argument('--cds', required=True, type=str, help='Bed table listing all CDS regions')
    parser.add_argument('--output', required=True, type=str, help='Output consensus')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    cds = pd.read_table(args.cds, header=None)
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        try:
            name
        except NameError:
            name = "empty_consensus"
        if name == "empty_consensus":
            print("WARNING: No reads available: there will be no information about Ns")
        else:
			cds_name = 
    