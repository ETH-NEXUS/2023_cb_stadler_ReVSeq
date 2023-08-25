######
## Script name: aggregate_table.py
## Date: August 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the assignment tables
##        for each sample in a plate and aggregates them in a plate-level
##        summary
######

import pandas as pd
import argparse, os, sys


def get_top_strain(inputdir, sample, dirname):
    sampledir = inputdir + "/" + sample
    if(dirname == "assign_virus"):
        strain = pd.read_table(sampledir + "/" + dirname + "/" + sample + "_substrain_count_table.tsv")
    elif(dirname == "validate_assignment"):
        strain = pd.read_table(sampledir + "/" + dirname + "/" + sample + "_validation.tsv")
    else:
        sys.exit("ERROR: unknown subdir to fetch the assignments")
    top_strain = strain.loc[strain['rpkm_proportions'] == max(strain['rpkm_proportions'])]
    if len(top_strain) > 1:
        if top_strain['coverage'][0] == 0:
            top_strain = top_strain.iloc[0]
            top_strain['name'] = ""
        else:
            sys.exit("Error: multiple highest strains with exactly the same rpkm proportion for sample: " + sample)
    return(top_strain)


def move_last_column_first(df):
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df


def create_line(df, linetype, sample):
    if linetype == "substrain":
        names = ["substrain_name", "aligned_reads_substrain", "reference_length", "rpkm_substrain", "rpkm_proportions_substrain", "outlier_substrain", "percentile_threshold_substrain"]
        df = df.iloc[:, 0:7]
    elif linetype == "strain":
        names = ["strain_name", "aligned_reads_strain", "rpkm_strain", "rpkm_proportions_strain", "outlier_strain", "percentile_threshold_strain", "coverage", "coverage_threshold", "qc_status", "panel_positive", "sample"]
    else:
        sys.exit("ERROR: invalid linetype")
    df.columns = names
    if linetype == "strain":
        df = move_last_column_first(df)
    return df


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inputdir', required=True, type=str, help='The main plate results directory from the ReVSeq pipeline')
    parser.add_argument('--assignment_subdir', required=True, type=str, help='The name assigned to the assignment subdirectory')
    parser.add_argument('--validation_subdir', required=True, type=str, help='The name assigned to the validation subdirectory')
    parser.add_argument('--qc_status', required=True, type=str, help='The status output text file from rule qualimap_filtered')
    parser.add_argument('--outfile', required=True, type=str, help='the path and filename for the output')

    args = parser.parse_args()
    ignorable_files = ['logs', 'multiqc', 'multiqc_filtered', 'package_results', 'complete.txt',  'aggregate']

    files = os.listdir(args.inputdir)
    samples = [ file for file in files if file not in ignorable_files]

    for sample in samples:
        top_substrain = get_top_strain(args.inputdir, sample, args.assignment_subdir)
        top_substrain = create_line(top_substrain, "substrain", sample)
        top_strain = get_top_strain(args.inputdir, sample, args.validation_subdir)
        top_strain = create_line(top_strain, "strain", sample)
        with open(args.inputdir + "/" + sample + "/" + args.qc_status + "/qc_status.txt") as f:
            qc_status = f.readline().strip()
        if qc_status == "PASS":
            try:
                aggregated_substrain.append(top_substrain)
            except NameError:
                aggregated_substrain = top_substrain
            try:
                aggregated_strain.append(top_strain)
            except NameError:
                aggregated_strain = top_strain
            aggregated_strain["qc_status"] = qc_status

    aggregated_table =  aggregated_substrain.merge(aggregated_strain, on="sample")
    aggregated_table.to_csv(args.outfile, sep="\t", float_format='%.5f')