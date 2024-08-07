######
## Script name: count_n_cds.py
## Date: March 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the 
##        the final consensus sequence and reports the status of the Ns
##        assigned to low-coverage regions.
##        The script returns the number of Ns per viral sequence and the fraction
##        of Ns for the reference sequence. This specific script is built to work
##        on CDS-based consensus sequences
######

import pandas as pd
import argparse, sys
from Bio import SeqIO


## Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='Count the number of Ns in a given consensus sequence', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, help='Input CDS-based consensus sequence')
    parser.add_argument('--bed', required=True, type=str, help='Bed reference table containing all annoated CDS regions')
    parser.add_argument('--output', required=True, type=str, help='Output table')
    parser.add_argument('--coverage', required=True, type=str, help='Bed reference table including coverage for all consensus positions')
    parser.add_argument('--threshold', required=True, type=int, help='Threshold used at consensus generation to mask low-coverage regions')
    parser.add_argument('--match_table', required=True, type=str, help='Matching table listing links between GenBank IDs and substrain names')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    bed = pd.read_table(args.bed, header=None, dtype={2:int, 3:int})
    cov = pd.read_table(args.coverage, header=None, dtype={1:int, 2:int, 3:int})
    matching_table = pd.read_table(args.match_table, header=None)
    
    info_df = pd.DataFrame(columns=["substrain", "CDS_name", "number_n", "reference_length", "fraction_n", "mean_cov_non_n_positions"])
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        tmp = name.split(":")[0]
        substrain = matching_table.loc[matching_table[0] == tmp][3].to_string(index=None).strip()
        try:
            name
        except NameError:
            name = "empty_consensus"
        if name == "empty_consensus":
            print("WARNING: No reads available: there will be no information about Ns")
        else:
            number_n = sequence.count("N")
            if len(sequence) == number_n:
                line = {"substrain": substrain, "CDS_name": fasta.description, "number_n": int(number_n), "reference_length": len(sequence), "fraction_n": 1, "mean_cov_non_n_positions": 0}
            else:
                start = int(name.split(":")[1].split("-")[0])
                end = int(name.split(":")[1].split("-")[1])
                fasta_cov = cov.loc[cov[0]==name.split(":")[0]]
                fasta_cov_start = (fasta_cov[1] <= start) & (fasta_cov[2] > start)
                if fasta_cov_start.sum() != 1:
                    sys.exit("ERROR: there are multiple locations in the coverage bed file matching the start of exon " + name)
                fasta_cov_end = (fasta_cov[1] < end) & (fasta_cov[2] >= end)
                if fasta_cov_end.sum() != 1:
                    sys.exit("ERROR: there are multiple locations in the coverage bed file matching the end of exon " + name)
                cov_start = fasta_cov_start.keys()[fasta_cov_start].to_list()
                cov_end = fasta_cov_end.keys()[fasta_cov_end].to_list()
                # If the data frame has length 1, iloc fails and we can just use fasta_cov as is
                if len(fasta_cov) != 1:
                    fasta_cov = cov.iloc[list(range(cov_start[0],cov_end[0]+1))]
                fasta_cov = fasta_cov.drop_duplicates()
                non_n_pos = [i for i, letter in enumerate(sequence) if (letter != "N")]
                full_cov = pd.DataFrame(columns={0,1,2})
                for index, row in fasta_cov.iterrows():
                    # we are not using a +1 in the range because the coverage is 0 based and goes 0-71 -> 71-73... With the + 1 we would have some duplicate positions
                    full_cov = full_cov.append(pd.DataFrame({0:row[0], 1:list(range(row[1], (row[2]))), 2:row[3]}))
                full_cov = full_cov.loc[full_cov[1] >= start]
                full_cov = full_cov.loc[full_cov[1] <= end]
                full_cov = full_cov.loc[full_cov[2] >= args.threshold ]
                if (len(full_cov) == 0):
                    avg = 0
                else:
                    length = len(full_cov)
                    avg = full_cov[2].sum() / length
                line = {"substrain": substrain, "CDS_name": fasta.description, "number_n": int(number_n), "reference_length": len(sequence), "fraction_n": 1, "mean_cov_non_n_positions": avg}
            info_df = info_df.append([line], ignore_index=True)
    if len(info_df) != 0:
        info_df["fraction_n"] = info_df["number_n"]/info_df["reference_length"]
    info_df = info_df.drop("reference_length", axis=1)
    info_df.to_csv(args.output, sep="\t")