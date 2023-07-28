######
## Script name: anonymize_sample_names.py
## Date: July 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq command that controls the revseq workflow.
##      Generates random, unique 6-characters IDs used to anonymize the sample
##      names provided by the sequencing lab.
##	    Replicates the raw data directory structure in a new directory using softlinks.
##      Creates the anonymization table with the matches between original IDs and anonymized IDs.
##      Creates the sample map file for the pipeline with the list of lanes available for each sample.
######


import pandas as pd
import sys, os
import argparse
import uuid
import base64
import unittest
import gzip

def truncated_uuid4():
    return str(uuid.uuid4())[:6]


def detect_empty_sample(samplename, anondir):
    # samples with 0 reads crash the pipeline and should be skipped
    # these samples will not be added to the samplemap to avoid crashing the procedure
    with gzip.open(anondir + "/" + samplename + "/" + samplename + "_L001_R1_001.fastq.gz", 'rb') as f:
        data = f.read(1)
    if len(data) == 0:
        print("Empty sample " + samplename + " will be excluded from the analysis")
        return [ samplename, "empty" ]
    else:
        return [ samplename, "not_empty" ]


def generate_new_ethids(anon_df, all_samples):
    if len(anon_df) == 0:
        print("Warning: empty anonymization table. This is just a notice, not an error.")
    newsamples = []
    for samplename in all_samples:
        if len(anon_df) > 0:
            if int(samplename.name) in anon_df["sample_name"].values:
                print("Sample already in anonymization table. Skipping.")
                continue
            else:
                newethid = truncated_uuid4()
                while (newethid in anon_df["ethid"]) and any(newethid == sublist[1] for sublist in newsamples):
                    newethid = truncated_uuid4()
                newsamples.append( [samplename.name, newethid] )
        else:
            newethid = truncated_uuid4()
            while any(newethid == sublist[1] for sublist in newsamples):
                newethid = truncated_uuid4()
            newsamples.append( [samplename.name, newethid] )
    newsamples = pd.DataFrame(newsamples)
    return newsamples


def link_anonimised_names(samplename, ethid, sampledir, anonymizeddir):
    print("Log: working on ethid: ", ethid)
    sampledir = sampledir + "/" + str(samplename)
    anonymizeddir = anonymizeddir + "/" + str(ethid)
    if not os.path.exists(sampledir):
        sys.exit("Error: there is no directory to be used as origin for the anonymized link for ethid " + str(ethid) + ". This should not happen. Please check if there are problems with the storage.\nPLEASE NOTE that all links for the previous ethids in this run have been already generated and need to be cleaned up before the next run. At this point, the anonymization matching table has not been updated yet")
    try:
        os.mkdir(anonymizeddir)
        for file in os.scandir(sampledir):
            temp_split = file.name.split("_")
            temp_split[0] = str(ethid)
            anonfile = "_".join([temp_split[0], temp_split[2], temp_split[3], temp_split[4]])
            os.symlink(sampledir + "/" + file.name, anonymizeddir + "/" + anonfile)
    except:
        for anonfile in os.scandir(anonymizeddir):
            os.remove(anonymizeddir + "/" + anonfile)
        os.rmdir(anonymizeddir)
        print("Warning: Failed to link the raw data to the anonymized folder for ethid " + str(ethid) + ". The sample will be skipped and not saved in the anonymization table")
        return [samplename, ethid, 1]
    return [samplename, ethid, 0]


def create_sample_map(samples, anondir):
    lanes = []
    for id in samples:
        rawfiles = os.scandir(anondir + "/" + id)
        lanes_partial = []
        for file in rawfiles:
            lanes_partial.append(file.name.split("_")[1])
        lanes_partial = list(set(lanes_partial))
        lanes_partial = [ [id, lane] for lane in lanes_partial]
        for lane in lanes_partial:
            lanes.append(lane)
    lanes = pd.DataFrame(lanes)
    lanes = lanes.rename(columns={0: "sample", 1: "lane"})
    return lanes


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--anontable', required=True, type=str, help='the table containing the full anonymization list')
    parser.add_argument('--sampledir', required=True, type=str, help='the directory containing all raw data samples, one folder per sample')
    parser.add_argument('--samplemapfile', required=True, type=str, help='the path and filename to the samplemap used by the pipeline')
    parser.add_argument('--anonymizeddir', required=True, type=str, help='the directory containing all anonymized samples, as links to the original data')
    parser.add_argument('--emptyfile', required=True, type=str, help='the path and filename to the file that contains all empty samples to report')

    args = parser.parse_args()
    
    anon = pd.read_table(args.anontable, header=0)
    all_samples = os.scandir(args.sampledir)

    newsamples = generate_new_ethids(anon, all_samples)

    linked = [ link_anonimised_names(sample, ethid, args.sampledir, args.anonymizeddir) for sample,ethid in zip(newsamples[0], newsamples[1]) ]

    processed = [ sample[0:2] for sample in linked if sample[2]==0 ]
    try:
        pd.DataFrame(processed).to_csv(args.anontable, sep="\t", mode="a", header=None, index=None)
    except:
        sys.exit("Error: could not write the anonymization table. All directories have been created and need to be deleted for the procedure to move forward\nTo know what directories to delete, check in the anonymized directory the ethids that are not included in the anonymised table and delete all of them")

    all_samples_status = [ detect_empty_sample(sample[1], args.anonymizeddir) for sample in processed ]
    empty_samples = [ sample[0] for sample in all_samples_status if sample[1] == "empty"]
    with open(args.emptyfile, 'w') as f:
        f.write("\n".join(empty_samples))

    not_empty_samples = [ sample[0] for sample in all_samples_status if sample[1] == "not_empty"]

    samplemap = create_sample_map(not_empty_samples, args.anonymizeddir)
    try:
        samplemap.to_csv(args.samplemapfile, sep="\t", mode="a", header=None, index=None)
    except:
        sys.exit("Error: could not write the samplemap table. All directories and new lines in the anonymization table have been created and need to be deleted for the procedure to move forward.")

    
