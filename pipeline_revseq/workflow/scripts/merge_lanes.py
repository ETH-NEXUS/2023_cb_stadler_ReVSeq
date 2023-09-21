######
## Script name: merge_lanes.py
## Date: August 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the pseudoanonymization table
##        and merges the multi-lane raw files while applying the pseudoanonymization
######

import argparse, shutil
import pandas as pd

# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--rawdir', required=True, type=str, help='The directory containing the raw data')
    parser.add_argument('--lanefiles', required=True, type=str, nargs='+', help='List containing all filenames for the sample under analysis')
    parser.add_argument('--pseudoanon_table', required=True, type=str, help='Location of the pseudoanonymization table')
    parser.add_argument('--outdir', required=True, type=str, help='the output directory')
    parser.add_argument('--outsuffix', required=True, type=str, help='Suffix to use, after the sample name, for the output file. It must include the file extension')

    args = parser.parse_args()

    matches = pd.read_table(args.pseudoanon_table)

    sample_name = int(args.lanefiles[0].split("/")[-1].split("_")[0])
    pseudoanon = matches.loc[matches['Sample number'] == sample_name]["ethid"]
    if len(pseudoanon) == 0:
        sys.exit("ERROR: cannot find a match for sample " + sample_name + " in the pseudoanonymization table.")
    if len(pseudoanon) > 1:
        sys.exit("ERROR: More than 1 match for sample " + sample_name + " in the pseudoanonymization table.")

    pseudoanon = pseudoanon.to_string(index=False)
    outfile = args.outdir + "/" + pseudoanon + args.outsuffix
    with open(outfile, 'wb') as wfp:
        for fn in args.lanefiles:
            with open(fn, 'rb') as rfp:
                shutil.copyfileobj(rfp, wfp)
