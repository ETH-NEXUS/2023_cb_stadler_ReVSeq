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
    top_strain = top_strain.drop(columns="Unnamed: 0", errors="ignore")
    if len(top_strain) > 1:
        if top_strain['aligned_reads'].sum() == 0:
            top_strain = top_strain.head(1)
            top_strain['name'] = ""
            top_strain['outlier'] = ""
            top_strain = top_strain.drop(columns="Unnamed: 0", errors="ignore")
        else:
            sys.exit("Error: multiple highest strains with exactly the same rpkm proportion for sample: " + sample)
    return(top_strain)


def get_alternative_tops(inputdir, sample, dirname, top_strain_name):
    if len(top_strain_name) == 0:
        return ""
    top_strain_name = top_strain_name[0]
    sampledir = inputdir + "/" + sample
    if(dirname == "assign_virus"):
        strain = pd.read_table(sampledir + "/" + dirname + "/" + sample + "_substrain_count_table.tsv")
        columns_to_use = ['name', 'rpkm_proportions', 'aligned_reads', "DP"]
    elif(dirname == "validate_assignment"):
        strain = pd.read_table(sampledir + "/" + dirname + "/" + sample + "_validation.tsv")
        columns_to_use = ['name', 'rpkm_proportions', 'aligned_reads']
    else:
        sys.exit("ERROR: unknown subdir to fetch the assignments")
    top_strain = strain.loc[strain['outlier'] == "*"]
    outliers = top_strain[columns_to_use].values.tolist()
    if len(outliers) == 1 or len(outliers) == 0:
        return ""
    matches = [ item for item in outliers if top_strain_name in item ]
    if len(matches) != 1:
        sys.exit("ERROR: multiple matches for the top strain. This should never happen")
    outliers.remove(matches[0])
    for item in outliers:
        if item == outliers[0]:
            if(dirname == "assign_virus"):
                prefix = " (DP:" + f'{item[3]:.5f}' + "; "
            elif(dirname == "validate_assignment"):
                prefix = " ("
            else:
                sys.exit("ERROR: unknown subdir to fetch the assignments")
            mystring = item[0] + prefix + "reads:" + f'{item[2]:.5f}' + "; rpkm_prop: " + f'{item[1]:.5f}' + ")"
        else:
            mystring = mystring + ", " + item[0] + " (DP:" + f'{item[3]:.5f}' + "; reads:" + f'{item[2]:.5f}' + " rpkm_prop: " + f'{item[1]:.5f}' + ")"
    return(mystring)


def create_line(df, linetype, sample, alternative_outliers):
    if linetype == "substrain":
        names = ["substrain_name", "aligned_reads_substrain", "reference_length", "DP", "rpkm_substrain", "rpkm_proportions_substrain", "outlier_substrain", "percentile_threshold_substrain", "readnum_threshold", "readnum_status", "DP_threshold", "DP_status", "alternative_outliers_substrain", "sample"]
        names_ordered = ["sample", "substrain_name", "aligned_reads_substrain", "reference_length", "DP", "rpkm_substrain", "rpkm_proportions_substrain", "outlier_substrain", "percentile_threshold_substrain", "alternative_outliers_substrain", "readnum_threshold", "readnum_status", "DP_threshold", "DP_status"]
        df['alternative_outliers_substrain'] = alternative_outliers
        df["sample"] = sample
    elif linetype == "strain":
        df['alternative_outliers_strain'] = alternative_outliers
        df["sample"] = sample
        names = ["strain_name", "aligned_reads_strain", "rpkm_strain", "rpkm_proportions_strain", "outlier_strain", "percentile_threshold_strain", "panel_positive", "alternative_outliers_strain", "sample"]
        names_ordered = ["sample", "strain_name", "aligned_reads_strain", "rpkm_strain", "rpkm_proportions_strain", "outlier_strain", "percentile_threshold_strain", "alternative_outliers_strain", "panel_positive"]
    else:
        sys.exit("ERROR: invalid linetype")
    df.columns = names
    df = df[names_ordered]
    return df


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inputdir', required=True, type=str, help='The main plate results directory from the ReVSeq pipeline')
    parser.add_argument('--assignment_subdir', required=True, type=str, help='The name assigned to the assignment subdirectory')
    parser.add_argument('--validation_subdir', required=True, type=str, help='The name assigned to the validation subdirectory')
    parser.add_argument('--outfile', required=True, type=str, help='the path and filename for the output')
    parser.add_argument('--sample_map', required=True, type=str, help='Sample map file containing all sample names')

    args = parser.parse_args()
    sample_map = pd.read_table(args.sample_map)
    samples = sample_map["sample"].unique()

    files = os.listdir(args.inputdir)
    samples = [file for file in files if file in samples]

    for sample in samples:
        top_substrain = get_top_strain(args.inputdir, sample, args.assignment_subdir)
        alternative_outliers_substrain = get_alternative_tops(args.inputdir, sample, args.assignment_subdir, top_substrain["name"].values)
        top_substrain = create_line(top_substrain, "substrain", sample, alternative_outliers_substrain)
        top_strain = get_top_strain(args.inputdir, sample, args.validation_subdir)
        alternative_outliers_strain = get_alternative_tops(args.inputdir, sample, args.validation_subdir, top_strain["name"].values[0])
        top_strain = create_line(top_strain, "strain", sample, alternative_outliers_strain)
        try:
            aggregated_substrain = pd.concat([aggregated_substrain, top_substrain], ignore_index=True)
        except NameError:
            aggregated_substrain = top_substrain
        try:
            aggregated_strain = pd.concat([aggregated_strain, top_strain], ignore_index=True)
        except NameError:
            aggregated_strain = top_strain

    aggregated_table =  aggregated_substrain.merge(aggregated_strain, on="sample")
    aggregated_table.to_csv(args.outfile, sep="\t", float_format='%.5f')