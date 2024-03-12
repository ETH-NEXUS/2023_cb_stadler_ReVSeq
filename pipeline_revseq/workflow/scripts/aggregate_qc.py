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
    top_strain = top_strain.drop(columns="Unnamed: 0", errors="ignore")
    if len(top_strain) > 1:
        if top_strain['coverage'][0] == 0:
            top_strain = top_strain.head(1)
            top_strain['name'] = ""
            top_strain['outlier'] = ""
            top_strain = top_strain.drop(columns="Unnamed: 0", errors="ignore")
        else:
            sys.exit("Error: multiple highest strains with exactly the same rpkm proportion for sample: " + sample)
    return(top_strain)


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='Create an aggregate table with all read counts detected through the pipeline steps', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inputdir', required=True, type=str, help='The main plate results directory from the ReVSeq pipeline')
    parser.add_argument('--bwa_subdir', required=True, type=str, help='The name assigned to the bwa rule subdirectory')
    parser.add_argument('--dehuman_subdir', required=True, type=str, help='The name assigned to the dehuman rule subdirectory')
    parser.add_argument('--filter_subdir', required=True, type=str, help='The name assigned to the filter_alignment rule subdirectory')
    parser.add_argument('--duplicates_subdir', required=True, type=str, help='The name assigned to the remove_duplicates subdirectory')
    parser.add_argument('--outfile', required=True, type=str, help='the path and filename for the output')
    parser.add_argument('--sample_map', required=True, type=str, help='Sample map file containing all sample names')


    args = parser.parse_args()
    sample_map = pd.read_table(args.sample_map)
    samples = sample_map["sample"].unique()

    files = os.listdir(args.inputdir)
    samples = [file for file in files if file in samples]

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
        if bwa_counts == '0':
            dehuman_percent = 0
        else:
            dehuman_percent = int(dehuman_counts)/int(bwa_counts)*100
        if dehuman_counts == '0':
            filter_percent = 0
        else:
            filter_percent = int(filter_counts)/int(dehuman_counts)*100
        if filter_counts == '0':
            duplicate_percent = 0
        else:
            duplicate_percent = int(duplicate_counts)/int(filter_counts)*100
        if bwa_counts_all == '0':
            duplicate_percent_all = 0
            bwa_percent = 0
        else:
            duplicate_percent_all = int(duplicate_counts)/int(bwa_counts_all)*100
            bwa_percent = int(bwa_counts)/int(bwa_counts_all)*100
        line = pd.DataFrame({"sample_name": sample, "all_reads": bwa_counts_all, "aligned_reads": bwa_counts, "aligned - % of all": bwa_percent, "dehumanized_reads": dehuman_counts, "dehumanized - % of aligned": dehuman_percent, "passed_filters_reads": filter_counts, "passed_filters - % of dehumanized": filter_percent, "deduplicated_reads": duplicate_counts, "deduplicated - % of filtered": duplicate_percent, "all_filters - % of total": duplicate_percent_all}, index=[0])
        try:
            aggregated_qc = aggregated_qc.append(line)
        except NameError:
            aggregated_qc = line

    aggregated_qc.to_csv(args.outfile, sep="\t", float_format='%.5f', index=None)