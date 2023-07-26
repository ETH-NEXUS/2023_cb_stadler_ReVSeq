from Bio import SeqIO, SeqUtils
import pandas as pd
import argparse

# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='plot coverage, diversity and output consensus sequence',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--refs', required=True, type=str, help='multifasta file with the viral reference sequences')
    parser.add_argument('--output', type=str, help='filename to use as output')

    refs = SeqIO.parse(args.refs, "fasta")
    gc_content = []
    for i in refs.records:
        gc_content.append([i.name, SeqUtils.GC(i)])
        gc_table = pd.DataFrame(gc_content)
        gc_table = gc_table.rename(columns={0: "reference", 1:"gc_content"})
        gc_table.to_csv(args.output, index=None)