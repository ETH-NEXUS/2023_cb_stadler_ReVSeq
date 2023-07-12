import glob, sys,os, argparse
from Bio import SeqIO, SeqRecord, Seq

# Script
if __name__ == '__main__':
    import pandas as pd

    # Parse input args
    parser = argparse.ArgumentParser(description='plot coverage, diversity and output consensus sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--refs', required=True, type=str, help='multifasta file with the viral reference sequences')
    parser.add_argument('--host_ref', required=True, type=str, help='fasta file with the host reference sequence')
    parser.add_argument('--output', type=str, help='filename to use as output')

    args = parser.parse_args()

    refs = SeqIO.parse(args.refs, "fasta")
    with open(args.output, "w") as output_handle:
        SeqIO.write(refs, output_handle, "fasta")

    with open(args.output, "a") as output_handle:                
        with open(args.host_ref) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                record.id = record.id + "_host"
                SeqIO.write(record, output_handle, "fasta")


