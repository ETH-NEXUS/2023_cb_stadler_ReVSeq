######
## Script name: aggregate_qc.py
## Date: August 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the QC read counts
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
    return(top_strain)


def move_last_column_first(df):
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    return df


def create_line(df, linetype, sample):
    if linetype == "substrain":
        names_dict = {"name": "substrain_name", "aligned_reads": "aligned_reads_substrain", "rpkm": "rpkm_substrain", "rpkm_proprtions": "rpkm_proportions_substrain", "outlier": "outlier_substrain", "percentile_threshold": "percentile_threshold_substrain"}
    elif linetype == "strain":
        names_dict = {"name": "strain_name", "aligned_reads": "aligned_reads_strain", "rpkm": "rpkm_strain", "rpkm_proprtions": "rpkm_proportions_strain", "outlier": "outlier_strain", "percentile_threshold": "percentile_threshold_strain"}
    else:
        sys.exit("ERROR: invalid linetype")
    df = df.rename(columns=names_dict)
    df['sample'] = sample
    df = move_last_column_first(df)
    return df


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inputdir', required=True, type=str, help='The main plate results directory from the ReVSeq pipeline')
    parser.add_argument('--bwa_subdir', required=True, type=str, help='The name assigned to the bwa rule subdirectory')
    parser.add_argument('--dehuman_subdir', required=True, type=str, help='The name assigned to the dehuman rule subdirectory')
    parser.add_argument('--filter_subdir', required=True, type=str, help='The name assigned to the filter_alignment rule subdirectory')
    parser.add_argument('--duplicates_subdir', required=True, type=str, help='The name assigned to the remove_duplicates subdirectory')
    parser.add_argument('--outfile', required=True, type=str, help='the path and filename for the output')

    args = parser.parse_args()
    ignorable_files = ['logs', 'multiqc', 'multiqc_filtered', 'package_results', 'complete.txt',  'aggregate']

    #add args here
    files = os.listdir(args.inputdir)
    samples = [ file for file in files if file not in ignorable_files]

    for sample in samples:
        sampledir = args.inputdir + "/" + sample
        with open(sampledir + "/" + args.bwa_subdir + "/" + sample + "_readcount_all.txt") as f:
            bwa_counts_all = f.readline().strip()
        with open(sampledir + "/" + args.bwa_subdir + "/" + sample + "_readcount.txt") as f:
            bwa_counts = f.readline().strip()
        with open(sampledir + "/" + args.filter_subdir + "/" + sample + "_readcount.txt") as f:
            filter_counts = f.readline().strip()
        with open(sampledir + "/" + args.dehuman_subdir + "/" + sample + "_readcount.txt") as f:
            dehuman_counts = f.readline().strip()
        with open(sampledir + "/" + args.duplicates_subdir + "/" + sample + "_readcount.txt") as f:
            duplicate_counts = f.readline().strip()
        line = pd.DataFrame({"sample_name": sample, "all_reads": bwa_counts_all, "aligned_reads": bwa_counts, "dehumanized_reads": dehuman_counts, "filtered_reads": filter_counts, "deduplicated_reads": duplicate_counts}, index=[0])
        try:
            aggregated_qc.append(line)
        except NameError:
            aggregated_qc = line

    aggregated_qc.to_csv(args.outfile, sep="\t", float_format='%.5f', index=None)