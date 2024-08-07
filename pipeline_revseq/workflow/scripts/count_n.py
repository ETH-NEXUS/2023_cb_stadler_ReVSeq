######
## Script name: count_n.py
## Date: March 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the 
##        the final consensus sequence and reports the status of the Ns
##        assigned to low-coverage regions.
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
    parser.add_argument('--coverage', required=True, type=str, help='Bed reference table including coverage for all consensus positions')
    parser.add_argument('--threshold', required=True, type=int, help='Threshold used at consensus generation to mask low-coverage regions')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    bed = pd.read_table(args.bed, header=None)
    cov = pd.read_table(args.coverage, header=None, dtype={1:int, 2:int, 3:int})

    info_df = pd.DataFrame(columns=["name", "number_n", "reference_length", "fraction_n", "mean_cov_non_n_positions"])
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        try:
            name
        except NameError:
            name = "empty_consensus"
        if name == "empty_consensus":
            print("WARNING: No reads available: there will be no information about Ns")
        else:
            fasta_cov = cov.loc[cov[0]==name]
            fasta_cov = fasta_cov.loc[fasta_cov[3] >= args.threshold]
            if (len(fasta_cov) == 0):
                avg = 0
            else:
                length = fasta_cov[2] - fasta_cov[1]
                avg = (length * fasta_cov[3]).sum() / length.sum()
            reference_length = bed[2][bed[0] == name].to_string(index=False).strip()
            name = bed[3][bed[0] == name].to_string(index=False).strip()
            number_n = sequence.count("N")
            line = {"name": name, "number_n": int(number_n), "reference_length": int(reference_length), "fraction_n": 1, "mean_cov_non_n_positions": avg}
            info_df = info_df.append([line], ignore_index=True)
    if len(info_df) != 0:
        info_df = info_df.groupby(["name"])[["number_n", "reference_length", "fraction_n", "mean_cov_non_n_positions"]].apply(lambda x : x.astype(float).sum())
        info_df["fraction_n"] = info_df["number_n"]/info_df["reference_length"]
    info_df = info_df.drop("reference_length", axis=1)
    info_df.to_csv(args.output, sep="\t")