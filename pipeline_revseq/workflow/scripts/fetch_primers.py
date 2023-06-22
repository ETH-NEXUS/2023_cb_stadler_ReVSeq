import glob, sys,os, argparse
import pandas as pd


def fetch_primers(primerfile, reference, output):
    import pandas as pd
    from Bio import SeqIO, Seq
    raw_primers = pd.read_csv(primerfile, sep='\t')
    ref = str(SeqIO.read(reference, 'fasta').seq)
    primers = {}
    for r, row in raw_primers.iterrows():
        start = ref.find(row.seq)
        if start<0:
            start = ref.find(Seq.reverse_complement(row.seq))
        if start>0:
            primers[row.name] = {"segment":segments[0], "name":row["name"], "seq":row.seq, "start":start, "end":start+len(row.seq)}
        else:
            print(f"row {row} failed")
    pd.DataFrame(primers).T.to_csv(output, sep='\t', index=False)

# Script
if __name__ == '__main__':
    import pandas as pd
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--primerfile', required=True, type=str, help='the file containing the list of primers')
    parser.add_argument('--reference', required=True, type=str, help='the reference file in FASTA format')
    parser.add_argument('--output', required=True, type=str, help='the name of the output file')

    args = parser.parse_args()
    stats = {}
    fetch_primers(args.primerfile, args.reference, args.output)
