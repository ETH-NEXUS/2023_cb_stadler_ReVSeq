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
    parser.add_argument('--bed', required=True, type=str, help='Bed reference table matching genbank accessions with scientific names')
    parser.add_argument('--output', required=True, type=str, help='Output table')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    bed = pd.read_table(args.bed, header=None)

    info_df = pd.DataFrame(columns=["name", "number_n", "reference_length", "fraction_n"])
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        try:
            name
        except NameError:
            name = "empty_consensus"
        if name == "empty_consensus":
            print("WARNING: No reads available: there will be no information about Ns")
        else:
            reference_length = bed[2][bed[0] == name].to_string(index=False).strip()
            name = bed[3][bed[0] == name].to_string(index=False).strip()
            number_n = sequence.count("N")
            line = {"name": name, "number_n": int(number_n), "reference_length": int(reference_length), "fraction_n": 1}
            info_df = info_df.append([line], ignore_index=True)
    if len(info_df) != 0:
        info_df = info_df.groupby(["name"])[["number_n", "reference_length"]].apply(lambda x : x.astype(int).sum())
        info_df["fraction_n"] = info_df["number_n"]/info_df["reference_length"]
    info_df = info_df.drop("reference_length", axis=1)
    info_df.to_csv(args.output, sep="\t")