######
## Script name: count_n.py
## Date: March 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the 
##        the final consensus sequence and reports the status of the Ns
##        assigned to low-coverage reasons.
##        The script returns the number of Ns per viral sequence and the fraction
##        of Ns for the reference sequence
######

import pandas as pd
import argparse
from Bio import SeqIO


## Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='Count the number of Ns in a given consensus sequence', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, help='Input consensus sequence')
    parser.add_argument('--output', required=True, type=str, help='Output table')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    with open(args.output, "a") as out_file:
        out_file.write("GenBank_accession\tnumber_n\tfraction_n")
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            number_n = sequence.count("N")
            fraction_n = number_n/len(sequence)
            out_file.write(name+"\t"+str(number_n)+"\t"+str(fraction_n))