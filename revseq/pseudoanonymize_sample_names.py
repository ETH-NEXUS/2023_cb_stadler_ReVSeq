######
## Script name: pseudoanonymize_sample_names.py
## Date: July 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq command that controls the revseq workflow.
##      Generates random, unique 6-characters IDs used to anonymize the sample
##      names provided by the sequencing lab.
##	    Replicates the raw data directory structure in a new directory using softlinks.
##      Creates the anonymization table with the matches between original IDs and pseudoanonymized IDs.
##      Creates the sample map file for the pipeline with the list of lanes available for each sample.
######


import pandas as pd
import sys, os, re
import argparse
import shortuuid
import base64
import unittest
import gzip

def generate_id(length=6):
    return shortuuid.ShortUUID().random(length=length)


def detect_empty_sample(samplename, outdir, plate):
    # samples with 0 reads crash the pipeline and should be skipped
    # these samples will not be added to the samplemap to avoid crashing the procedure
    inputdir = outdir + "/" + plate + "/" + samplename 
    file = os.path.join(inputdir, os.listdir(inputdir)[0])
    with gzip.open(file, 'rb') as f:
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
        newethid = generate_id()
        if len(anon_df) > 0:
            while (newethid in anon_df) and any(newethid == sublist[1] for sublist in newsamples):
                newethid = generate_id()
            newsamples.append( [samplename, newethid] )
        else:
            while any(newethid == sublist[1] for sublist in newsamples):
                newethid = generate_id()
            newsamples.append( [samplename, newethid] )
    newsamples = pd.DataFrame(newsamples)
    return newsamples


def link_pseudoanonimised_names(samplename, ethid, sampledir, anonymizeddir, plate):
    print("Log: working on ethid: ", ethid)
    outdir = anonymizeddir + "/" + plate + "/" + ethid
    try:
        os.mkdir(outdir)
    except:
        sys.exit("ERROR: cannot create the pseudo-anonymized directory")
    fastqs = os.listdir(sampledir)
    regex = r"\A" + re.escape(str(samplename)) + r".+R.+.fastq.gz\Z"
    to_link = [rawfile for rawfile in fastqs if any(re.findall(regex, rawfile))]
    for rawfile in to_link:
        outfile = rawfile.replace(str(samplename), ethid)
        outfile = re.sub(r"_S\d{1,2}_", "_", outfile)
        try:
            os.symlink(args.mirrordir + "/" + str(rawfile), outdir + "/" + outfile)
        except:
            sys.exit("ERROR: cannot create the symlink for sample " + ethid)


def create_plate_directory(platename, outdir):
    print("Creating directory for plate: " + platename)
    try:
        os.mkdir(outdir + "/" + platename)
    except FileExistsError:
        sys.exit("ERROR: the folder " + outdir + "/" + platename, " already exisits")
    return outdir + "/" + platename


def create_sample_directories(platename, outdir, samples):
    print("Creating sample directories for plate: " + platename)
    for sample in samples:
        try:
            os.mkdir(outdir + "/" + platename + "/" + sample)
        except FileExistsError:
            sys.exit("ERROR: the folder " + outdir + "/" + platename, + "/" + sample + " already exisits")
    return


def create_sample_map(samples, outdir):
    lanes = []
    for id in samples:
        rawfiles = os.scandir(outdir + "/" + plate + "/" + id)
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

#def verify_if_new_plate(metadatadir, processed_plates):
#    metadata_subdirs = os.listdir(metadatadir)
#    metadata_files = [ metadatadir + "/" + subdir + "/" + os.listdir(metadatadir + "/" + subdir)[0] for subdir in metadata_subdirs ]
#    new_plates = []
#    for file in metadata_files:
#        with open(file) as f:
#            metadata = f.read().splitlines()
#        try:
#            plate = metadata[1].split(";")[4]
#        except indexError:
#            sys.exit("ERROR: the metadata table " + file + " is either empty or truncated")
#        if plate not in processed_plates:
#            new_plates.append([plate, file])
#    return new_plates


def verify_if_new_plate(mirrordir, mirrored_plates, processed_plates):
    new_plates_names = [ plate for plate in mirrored_plates if plate not in processed_plates ]
    new_plates = []
    for plate in new_plates_names:
        platedir = mirrordir + "/" + plate
        metadata_file = [ file for file in os.listdir(platedir) if ".csv" in file]
        metadata_file = [ platedir + "/" + file for file in metadata_file if "seq_" in file ]
        if len(metadata_file) != 1:
            sys.exit("ERROR: found " + str(len(metadata_file)) + " metadata files in mirrored directory " + platedir)
        metadata_file = metadata_file[0]
        with open(metadata_file) as f:
            metadata = f.read().splitlines()
        try:
            metadata_plate = metadata[1].split(";")[4]
            metadata_plate = metadata_plate.replace('"', '')
        except indexError:
            sys.exit("ERROR: the metadata table " + file + " is either empty or truncated")
        if metadata_plate != plate:
            sys.exit("ERROR: the metadata table " + file + " lists a different plate than the directory name")
        new_plates.append([plate,metadata_file])
    return new_plates


def verify_if_complete_plate(new_plates, mirrordir):
    completeness = {}
    # find how many samples we have in the mirror for the new plates
    for plate,file in new_plates:
        mirrorsamples = os.listdir(mirrordir + "/" + plate)
        samples = pd.read_csv(file, sep=";")['Sample number'].to_list()
        r1 = []
        r2 = []
        for sample in samples:
            regexr1 = r"\A" + re.escape(str(sample)) + r".+R1.+.fastq.gz\Z"
            regexr2 = r"\A" + re.escape(str(sample)) + r".+R2.+.fastq.gz\Z"
            r1.extend(set([sample for fastq in mirrorsamples if any(re.findall(regexr1, fastq))]))
            r2.extend(set([sample for fastq in mirrorsamples if any(re.findall(regexr2, fastq))]))
        if len(r1) == len(samples) and len(r2) == len(samples):
            completeness[plate] = [ r1, file ]
        else:
            completeness[plate] = None
    complete = {key:val for key, val in completeness.items() if val is not None}
    return complete


def get_pseudoanon_names(outdir):
    plates = [ outdir + "/" + plate for plate in os.listdir(outdir) ]
    anon = []
    for plate in plates:
        anon.extend(os.listdir(plate))
    return anon


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outdir', required=True, type=str, help='the path to the output directory')
    #parser.add_argument('--metadatadir', required=True, type=str, help='the path to the directory where the metadata are mirrored')
    parser.add_argument('--mirrordir', required=True, type=str, help='the directory containing all raw data samples, one folder per sample')

    args = parser.parse_args()

    processed_plates = os.listdir(args.outdir)
    if len(processed_plates) == 0:
        print("Warning: no processed plates found. If this is not the very first run of the workflow please check the configuration")
    
    mirrored_plates = os.listdir(args.mirrordir)
    if len(mirrored_plates) == 0:
        print("Warning: no mirrored plates found. Please check the configuration")
    
    anon = get_pseudoanon_names(args.outdir)

    #new_plates = verify_if_new_plate(args.metadatadir, processed_plates)
    new_plates = verify_if_new_plate(args.mirrordir, mirrored_plates, processed_plates)
    if len(new_plates) == 0:
        print("No new plate found.")
        print("NEW = the plate has not been pseudo-anonymized before")
        sys.exit()

    new_complete_plates = verify_if_complete_plate(new_plates, args.mirrordir)

    if len(new_complete_plates) == 0:
        print("There are new plates but none appears to be complete yet")
        print("COMPLETE = the raw data in the mirror show at least one R1 and one R2 file for each sample listed in the metadata.")
        sys.exit()
    
    print("Found new plates to process: " + str(new_plates))
    print("Beginning pseudo-anonymization")

    for plate,data in new_complete_plates.items():
        sample = data[0]
        file = data[1]
        newsamples = generate_new_ethids(anon, sample)
        newsamples = newsamples.rename(columns={0: "Sample number", 1: "ethid"})
        if len(newsamples.index) > 0:
            create_plate_directory(plate, args.outdir)
            create_sample_directories(plate, args.outdir, newsamples["ethid"].to_list())
            try:
                pd.DataFrame(newsamples).to_csv(outdir + "/" + plate +  "/" + plate + "_pseudoanon_table.tsv", sep="\t", index=None)
            except:
                sys.exit("Error: could not write the pseudoanonymization table. All directories have been created and need to be deleted for the procedure to move forward\nTo know what directories to delete, check in " + args.anonymizeddir + " the ethids that are not present in the pseudoanonymised table and delete all of them")
            all_samples_status = [ detect_empty_sample(sample, args.outdir, plate) for sample in newsamples["ethid"].to_list() ]
            empty_samples = [ sample[0] for sample in all_samples_status if sample[1] == "empty"]
            with open(outdir + "/" + plate + "_empty_samples.txt", 'w') as f:
                f.write("\n".join(empty_samples))

            metadata = pd.read_csv(file, sep=";")
            metadata = pd.merge(metadata, newsamples, on="Sample number")
            metadata = metadata.rename(columns={"Spital: Ambulant/Stanionär": "treatment_type"})
            metadata.to_csv(outdir + "/" + plate + "_metadata.csv", sep=";", index=None)

            not_empty_samples = [ sample[0] for sample in all_samples_status if sample[1] == "not_empty"]

            samplemap = create_sample_map(not_empty_samples, args.outdir)
            samplemap.rename(columns={0: "sample", 1: "lane"})
            try:
                samplemap.to_csv(outdir + "/" + plate + "_samplemap.tsv", sep="\t", index=None)
            except:
                sys.exit("Error: could not write the samplemap table. All directories and new lines in the pseudoanonymization table have been created and need to be deleted for the procedure to move forward.")
        else:
            print("Warning: all new samples appear to be already in the pseudo-anonymization table.")
            print("This is usually caused by samples that have been already pseudo-anonymized, but never processed by the pipeline.")
            print ("Skipping pseudo-anonymization")
    
