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
import argparse, os, sys, math


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


def get_alternative_tops(inputdir, sample, dirname, top_strain_name, filter_alternatives, readnum_threshold, coverage_threshold):
    if len(top_strain_name) == 0:
        return ""
    top_strain_name = top_strain_name[0]
    sampledir = inputdir + "/" + sample
    if(dirname == "assign_virus"):
        strain = pd.read_table(sampledir + "/" + dirname + "/" + sample + "_substrain_count_table.tsv")
        columns_to_use = ['name', 'rpkm_proportions', 'aligned_reads', "coverage"]
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
        if filter_alternatives:
            if (item[1] < readnum_threshold) or (item[3] < coverage_threshold):
                print("Filtered out alternative outlier "+str(item)+" for not hitting the read number or coverage thresholds")
                continue
        try:
            mystring = mystring + ", " + item[0] + " (cov:" + f'{item[3]:.5f}' + "; reads:" + f'{item[2]:.0f}' + " rpkm_prop: " + f'{item[1]:.5f}' + ")"
        except NameError:
            if(dirname == "assign_virus"):
                prefix = " (cov:" + f'{item[3]:.5f}' + "; "
            elif(dirname == "validate_assignment"):
                prefix = " ("
            else:
                sys.exit("ERROR: unknown subdir to fetch the assignments")
            mystring = item[0] + prefix + "reads:" + f'{item[2]:.0f}' + "; rpkm_prop: " + f'{item[1]:.5f}' + ")"
    try:
        mystring
    except NameError:
        return ""
    return(mystring)


def create_line(df, linetype, sample, alternative_outliers):
    if linetype == "substrain":
        names = ["substrain_name", "aligned_reads_substrain", "reference_length", "coverage", "rpkm_substrain", "rpkm_proportions_substrain", "outlier_substrain", "percentile_threshold_substrain", "readnum_threshold", "readnum_status", "coverage_threshold", "coverage_status", "taxon_id", "scientific_name", "DP1", "DP2", "DP20", "alternative_outliers_substrain", "sample"]
        names_ordered = ["sample", "substrain_name", "taxon_id", "scientific_name", "reference_length", "aligned_reads_substrain", "readnum_threshold", "readnum_status", "coverage", "coverage_threshold", "coverage_status", "DP1", "DP2", "DP20", "rpkm_substrain", "rpkm_proportions_substrain", "outlier_substrain", "percentile_threshold_substrain", "alternative_outliers_substrain"]
        df['alternative_outliers_substrain'] = alternative_outliers
        df["sample"] = sample
    elif linetype == "strain":
        df['alternative_outliers_strain'] = alternative_outliers
        df["sample"] = sample
        names = ["strain_name", "aligned_reads_strain", "rpkm_strain", "rpkm_proportions_strain", "outlier_strain", "percentile_threshold_strain", "panel_match",  "panel_assignment", "alternative_outliers_strain", "sample"]
        names_ordered = ["sample", "strain_name", "aligned_reads_strain", "rpkm_strain", "rpkm_proportions_strain", "outlier_strain", "percentile_threshold_strain", "alternative_outliers_strain", "panel_match", "panel_assignment"]
    else:
        sys.exit("ERROR: invalid linetype")
    df.columns = names
    df = df[names_ordered]
    return df


def get_assigned_panel(metadata, top_strain, match_table):
    metadata = metadata.loc[metadata['ethid'] == sample ]
    metadata = metadata.to_dict(orient='list')
    positive = []
    for key,value in metadata.items():
        value = value[0]
        try:
            value = float(value)
        except (ValueError, KeyError) as error:
            continue	
        # Current metadata has an empty cell, nan or "deleted" for missing values; 0 or -1 for negatives; 4 or a positive integer >4 for positives
        if (math.isnan(value) or (value <= 0) or (key == 'Sample number') or (key == "Aufnahmenummer")):
            continue
        if value > 0:
            positive.append(key.split(" ")[0])
    positive_name = []
    for virus in positive:
        positive_name.append(match_table.loc[match_table["panel_name"] == virus]["strain_name"].to_string(index=False))
    return positive_name


def test_qc(s):
    if (s["readnum_status"] == "SUCCESS") and (s["coverage_status"] == "SUCCESS"):
        if (s["alternative_outliers_substrain"] == "") and (s["alternative_outliers_strain"] == ""):
            return "Quality OK, virus assignment unique"
        elif (s["alternative_outliers_substrain"] != "") and (s["alternative_outliers_strain"] == ""):
            return "Quality OK, substrain ambiguous, strain unique"
        elif (s["alternative_outliers_substrain"] == "") and (s["alternative_outliers_strain"] != ""):
            return "Quality OK, substrain unique, strain ambiguous"
        elif (s["alternative_outliers_substrain"] != "") and (s["alternative_outliers_strain"] != ""):
            return "Quality OK, substrain and strain ambiguous"
    if (s["readnum_status"] == "FAILED") and (s["coverage_status"] == "SUCCESS"):
        return "Not enough evidence"
    if (s["readnum_status"] == "SUCCESS") and (s["coverage_status"] == "FAILED"):
        return "not enough evidence"
    if (s["readnum_status"] == "FAILED") and (s["coverage_status"] == "FAILED"):
        return "quality fail"


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='Aggregate the virus assignments in a single summary table', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inputdir', required=True, type=str, help='The main plate results directory from the ReVSeq pipeline')
    parser.add_argument('--assignment_subdir', required=True, type=str, help='The name assigned to the assignment subdirectory')
    parser.add_argument('--validation_subdir', required=True, type=str, help='The name assigned to the validation subdirectory')
    parser.add_argument('--outfile', required=True, type=str, help='the path and filename for the output')
    parser.add_argument('--sample_map', required=True, type=str, help='Sample map file containing all sample names')
    parser.add_argument('--pseudoanon_metadata', required=True, type=str, help='Metadata table including the pseudonymized sample names')
    parser.add_argument('--match_table', required=True, type=str, help='TSV matching the panel code with the common virus name used in the pipeline')
    parser.add_argument('--filter_alternatives', required=False, type=bool, default=False, help='Should we filter the alternative outliers using the same readnum and coverage thresholds used for the QC of the top outlier?')
    parser.add_argument('--readnum_threshold', required=True, type=int, help='The read number threshold to apply to the alternative outliers')
    parser.add_argument('--coverage_threshold', required=True, type=int, help='The coverage threshold to apply to the alternative outliers')
    parser.add_argument('--controls', required=True, type=str, help='The strings necessary to detect the controls. Comma separated')

    args = parser.parse_args()
    sample_map = pd.read_table(args.sample_map)
    samples = sample_map["sample"].unique()

    files = os.listdir(args.inputdir)
    samples = [file for file in files if file in samples]

    metadata = pd.read_csv(args.pseudoanon_metadata,sep=";")
    
    match_table = pd.read_csv(args.match_table,sep=",", header=0)

    ctrl_strings = args.controls.split(",")

    for sample in samples:
        top_substrain = get_top_strain(args.inputdir, sample, args.assignment_subdir)
        alternative_outliers_substrain = get_alternative_tops(args.inputdir, sample, args.assignment_subdir, top_substrain["name"].values, args.filter_alternatives, args.readnum_threshold, args.coverage_threshold)
        top_substrain = create_line(top_substrain, "substrain", sample, alternative_outliers_substrain)
        top_strain = get_top_strain(args.inputdir, sample, args.validation_subdir)
        alternative_outliers_strain = get_alternative_tops(args.inputdir, sample, args.validation_subdir, top_strain["name"].values[0], args.filter_alternatives, args.readnum_threshold, args.coverage_threshold)
        ctrl_type = ""
        for string in ctrl_strings:
            if string in sample:
                ctrl_type = string
        print(ctrl_type)
        if ctrl_type == "":
            top_strain["panel_assignment"] = ';'.join(get_assigned_panel(metadata, top_strain, match_table))
        else:
            top_strain["panel_assignment"] = ";".join("")
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
    aggregated_table["qc_status"] = aggregated_table.apply(test_qc, axis=1)
    aggregated_table.to_csv(args.outfile, sep="\t", float_format='%.5f')