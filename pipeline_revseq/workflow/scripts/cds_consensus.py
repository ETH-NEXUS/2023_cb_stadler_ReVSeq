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
import argparse, sys
from Bio import SeqIO
from Bio.Seq import Seq


## Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='Count the number of Ns in a given consensus sequence', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', required=True, type=str, help='Input consensus sequence')
    parser.add_argument('--cds', required=True, type=str, help='Bed table listing all CDS regions')
    parser.add_argument('--output', required=True, type=str, help='Output consensus')

    args = parser.parse_args()

    fasta_sequences = SeqIO.parse(open(args.input),'fasta')
    cds = pd.read_table(args.cds, header=None, dtype={2:int, 3:int})
    all_fasta_cds = []
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        try:
            name
        except NameError:
            name = "empty_consensus"
        if name == "empty_consensus":
            print("WARNING: No reads available: there will be no information about Ns")
        else:
            current_fasta_cds = cds.loc[cds[0]==name]
            if len(current_fasta_cds) == 0:
                current_fasta_cds = cds.loc[cds[0]==name+".1"]
                if len(current_fasta_cds) == 0:
                    print("ERROR: Cannot find CDS information for virus ID " + name)
            for index, row in current_fasta_cds.iterrows():
                if row[2] > len(fasta):
                    print("ERROR: CDS " + row.to_string() + " ends later than the available consensus length of " + str(len(fasta)) + " for virus ID " + name)
                current_fasta = fasta[row[1]:(row[2]+1)]
                current_fasta.id = current_fasta.id + " CDS " + str(row[1]) + "-" + str(row[2])
                all_fasta_cds.append(current_fasta)
    with open(args.output, "w") as output_handle:
        SeqIO.write(all_fasta_cds, output_handle, "fasta")
    