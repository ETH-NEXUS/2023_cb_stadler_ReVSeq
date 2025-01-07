######
## Script name: chr_file.py
## Date: December 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the assignment tables
##        for each sample in a plate, retrieves the top assignment and
##        create the chromosome file necessary for the uploads
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
    strain=strain.rename(columns={"Unnamed: 0": "name"})
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


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='Aggregate the virus assignments in a single summary table', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--inputdir', required=True, type=str, help='The main plate results directory from the ReVSeq pipeline')
    parser.add_argument('--assignment_subdir', required=True, type=str, help='The name assigned to the assignment subdirectory')
    parser.add_argument('--sample', required=True, type=str, help='The sample name')
    parser.add_argument('--k2_bed', required=True, type=str, help='Bed file containing all the viruses detected by Kraken2 plus the panel base. Created by pipeline step merge_regions')
    parser.add_argument('--out', required=True, type=str, help='The output file name')
    
args = parser.parse_args()

files = os.listdir(args.inputdir)
k2bed = pd.read_table(args.k2_bed, header=None, names=['id','start', 'end', 'name'])

top_substrain = get_top_strain(args.inputdir, args.sample, args.assignment_subdir)
if len(top_substrain) == 0:
    print("No top strain found for sample " + args.sample)
    with open(args.out, "w") as file:
        line = "\t\t"
        file.write(line)
else:
    name = top_substrain.iloc[0]['name']
    if "Influenza" in name:
        ids = k2bed.loc[k2bed['name'] == top_substrain.iloc[0]['name'], 'id']
        with open(args.out, "a") as file:
            for id in ids:
                line = id + "\t" + id + "\t" + "segmented\n"
                file.write(line)
    else:
        id = k2bed.loc[k2bed['name'] == top_substrain.iloc[0]['name'], 'id'].item()
        with open(args.out, "w") as file:
            line = id + "\t" + id + "\t" + "chromosome"
            file.write(line)

