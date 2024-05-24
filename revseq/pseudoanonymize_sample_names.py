######
## Script name: pseudoanonymize_sample_names.py
## Date: July 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq command that controls the revseq workflow.
##      Generates random, unique 6-characters IDs used to anonymize the sample
##      names provided by the sequencing lab.
##	    Replicates the raw data directory structure in a new directory using softlinks.
##      Creates the anonymization table with the matches between original IDs and pseudonymized IDs.
##      Creates the sample map file for the pipeline with the list of lanes available for each sample.
######


import pandas as pd
import sys, os, re, argparse, shortuuid, base64, unittest, gzip, glob

ctrls_error_msg = "ERROR: unexpected control file format. Expected format [KOpos|KOneg]_S[0-9][0-9]_L0[0-9][0-9]_R[1|2]_001.fastq.gz"

def generate_id(length=6):
    return shortuuid.ShortUUID().random(length=length)


def detect_empty_sample(samplename, mirrordir, plate):
    # samples with 0 reads crash the pipeline and should be skipped
    # these samples will not be added to the samplemap to avoid crashing the procedure
    samplename = str(samplename)
    inputdir = mirrordir + "/" + plate + "/"
    #file = os.path.join(inputdir, os.listdir(inputdir)[0])
    file = [ file for file in os.listdir(inputdir) if samplename in file ]
    if len(file) == 0:
        sys.exit("ERROR: cannot find the input fastq files in the mirror for sample " + samplename)
    file = file[0]
    with gzip.open(inputdir + "/" + file, 'rb') as f:
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
        if len(str(samplename).split("-")) == 3:
            newethid = "MouthWash-"+newethid
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


def find_controls(ctrl_neg_prefix, ctrl_pos_prefix, mirrordir, plate):
    ctrl_pos = [ glob.glob(mirrordir+"/"+plate+"/"+pos+"*") for pos in ctrl_pos_prefix.split(",") ]
    ctrl_pos = [ x for xs in ctrl_pos for x in xs ]
    ctrl_neg = [ glob.glob(mirrordir+"/"+plate+"/"+neg+"*") for neg in ctrl_neg_prefix.split(",") ]
    ctrl_neg = [ x for xs in ctrl_neg for x in xs ]
    all_ctrls = []
    all_ctrls.extend(ctrl_pos)
    all_ctrls.extend(ctrl_neg)
    ctrl_count = []
    for file in all_ctrls:
        filename = file.split("/")[-1]
        filename_pieces = filename.split("_")
        if ("L0" not in filename_pieces[2]) or (("R1" not in filename_pieces[3]) and ("R2" not in filename_pieces[3])):
            sys.exit(ctrls_error_msg)
        ctrl_count.append(filename_pieces[0])
    ctrl_count = set(ctrl_count)
    print("Found "+str(len((ctrl_count)))+" controls: "+str(ctrl_count))
    ctrl_name = [ element for element in ctrl_count]
    ctrl_ethid = [ plate+"-"+element.replace("_", "-") for element in ctrl_count ]
    ctrls = pd.DataFrame({"Sample number": ctrl_name, "ethid": ctrl_ethid})
    return ctrls


def link_pseudoanonimised_names(samplename, ethid, sampledir, anonymizeddir, plate):
    print("Log: working on ethid: ", ethid)
    outdir = anonymizeddir + "/" + plate + "/anonymized/"
    if not os.path.isdir(outdir):
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
            os.symlink(sampledir + "/" + str(rawfile), outdir + "/" + outfile)
        except:
            sys.exit("ERROR: cannot create the symlink for sample " + ethid)


def create_plate_directory(platename, outdir):
    print("Creating directory for plate: " + platename)
    try:
        os.mkdir(outdir + "/" + platename)
    except FileExistsError:
        sys.exit("ERROR: the folder " + outdir + "/" + platename, " already exisits")
    return outdir + "/" + platename


def create_sample_directories(platename, outdir, pseudoanon_table, mirrordir):
    print("Creating sample directories for plate: " + platename)
    for sample in pseudoanon_table["Sample number"].to_list():
        try:
            os.mkdir(outdir + "/" + platename + "/" + sample)
        except FileExistsError:
            sys.exit("ERROR: the folder " + outdir + "/" + platename, + "/" + sample + " already exisits")
    return


def create_sample_map(samples, outdir):
    lanes = []
    for id in samples:
        lanes_partial = []
        rawfiles = os.scandir(outdir + "/" + plate + "/anonymized")
        for file in rawfiles:
            lanes_partial.append(file.name.split("_")[1])
        lanes_partial = list(set(lanes_partial))
        lanes_partial = [ [id, lane] for lane in lanes_partial]
        for lane in lanes_partial:
            lanes.append(lane)
    lanes = pd.DataFrame(lanes)
    lanes = lanes.rename(columns={0: "sample", 1: "lane"})
    return lanes


def load_metadata(mirrordir, metadata_string, plate):
    platedir = mirrordir + "/" + plate
    metadata_file = [ file for file in os.listdir(platedir) if ".csv" in file]
    metadata_file = [ platedir + "/" + file for file in metadata_file if metadata_string in file ]
    if len(metadata_file) > 1:
        sys.exit("ERROR: found " + str(len(metadata_file)) + " metadata files in mirrored directory " + platedir)
    if len(metadata_file) == 0:
        print("WARNING: no metadata found for plate " + plate + ". Skipping...")
        return 0
    metadata_file = metadata_file[0]
    with open(metadata_file) as f:
        metadata = f.read().splitlines()
    return [metadata, metadata_file]


def verify_if_new_plate(mirrordir, mirrored_plates, processed_plates, metadata_string):
    new_plates_names = [ plate for plate in mirrored_plates if plate not in processed_plates ]
    new_plates = []
    for plate in new_plates_names:
        metadata = load_metadata(mirrordir, metadata_string, plate)
        if metadata == 0:
            continue
        metadata_file = metadata[1]
        metadata = metadata[0]
        try:
            metadata_plate = metadata[1].split(";")[4]
            metadata_plate = metadata_plate.replace('"', '')
        except indexError:
            sys.exit("ERROR: the metadata table for plate " + plate + " is either empty or truncated")
        if metadata_plate != plate:
            sys.exit("ERROR: the metadata table for plate " + plate + " lists a different plate than the directory name")
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


def find_ethid(sample_number, match_table):
    return match_table.loc[match_table["Sample number"]==str(sample_number)]["ethid"].to_string(index=False).strip()


def detect_mouthwash_samples(new_complete_plates, mirrordir):
    complete_plates = {}
    for plate,data in new_complete_plates.items():
        files = os.listdir(mirrordir + "/" + plate)
        for sample in files:
            substrings = sample.split("-")
            if len(substrings) != 3 or len(substrings[0]) != 1 or substrings[0].isnumeric() or len(substrings[1]) != 1 or (not substrings[1].isnumeric()):
                continue
            substrings[2] = substrings[2].split("_")[0]
            name = "-".join(substrings)
            if name not in data[0]:
                data[0].append(name)
        complete_plates[plate] = data
    return complete_plates


# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--outdir', required=True, type=str, help='the path to the output directory')
    #parser.add_argument('--metadatalist', required=True, type=str, help='Tab-delimited table containing three columns: "metadata_string" containing the string to use to fetch the metadata files, "id_column" containing the name of the column containing the IDs, "sep" containing the separator used in the file')
    parser.add_argument('--mirrordir', required=True, type=str, help='the directory containing all raw data samples, one folder per sample')
    parser.add_argument('--ctrl_neg_prefix', required=True, type=str, help='the prefixes to use to detect negative controls, comma separated')
    parser.add_argument('--ctrl_pos_prefix', required=True, type=str, help='the prefixes to use to detect positive controls, comma separated')
    parser.add_argument('--metadata_string', required=True, type=str, help='the string to use to find the metadata files')

    args = parser.parse_args()

    processed_plates = os.listdir(args.outdir)
    if len(processed_plates) == 0:
        print("Warning: no processed plates found. If this is not the very first run of the workflow please check the configuration")
    
    mirrored_plates = os.listdir(args.mirrordir)
    if len(mirrored_plates) == 0:
        print("Warning: no mirrored plates found. Please check the configuration")
    
    anon = get_pseudoanon_names(args.outdir)

    #new_plates = verify_if_new_plate(args.metadatadir, processed_plates)
    new_plates = verify_if_new_plate(args.mirrordir, mirrored_plates, processed_plates, args.metadata_string)
    if len(new_plates) == 0:
        print("No new plate found.")
        print("NEW = the plate has not been pseudo-anonymized before")
        sys.exit()

    new_complete_plates = verify_if_complete_plate(new_plates, args.mirrordir)
    complete_plates = detect_mouthwash_samples(new_complete_plates, args.mirrordir)

    if len(complete_plates) == 0:
        print("There are new plates but none appears to be complete yet")
        print("COMPLETE = the raw data in the mirror show at least one R1 and one R2 file for each sample listed in the metadata.")
        sys.exit()
    
    print("Found new plates to process: " + str(complete_plates.keys()))
    print("Beginning pseudo-anonymization")

    for plate,data in complete_plates.items():
        sample = data[0]
        file = data[1]
        newsamples = generate_new_ethids(anon, sample)
        newsamples = newsamples.rename(columns={0: "Sample number", 1: "ethid"})
        newsamples["Sample number"] = newsamples["Sample number"].astype(str)
        newsamples = pd.concat([newsamples, find_controls(args.ctrl_neg_prefix, args.ctrl_pos_prefix, args.mirrordir, plate)])
        sample = newsamples["Sample number"].tolist()
        if len(newsamples.index) > 0:
            create_plate_directory(plate, args.outdir)
            for s in sample:
                link_pseudoanonimised_names(s, find_ethid(s, newsamples), args.mirrordir+"/"+plate, args.outdir, plate)
            try:
                pd.DataFrame(newsamples).to_csv(args.outdir + "/" + plate +  "/" + plate + "_pseudoanon_table.tsv", sep="\t", index=None)
            except:
                sys.exit("Error: could not write the pseudoanonymization table. All directories have been created and need to be deleted for the procedure to move forward\nTo know what directories to delete, check in " + args.outdir + " the ethids that are not present in the pseudoanonymised table and delete all of them")
            all_samples_status = [ detect_empty_sample(sample, args.mirrordir, plate) for sample in newsamples["Sample number"].to_list() ]
            empty_samples = [ find_ethid(sample[0], newsamples) for sample in all_samples_status if sample[1] == "empty"]
            with open(args.outdir + "/" + plate + "/" + plate + "_empty_samples.txt", 'w') as f:
                f.write("\n".join(empty_samples))

            metadata = pd.read_csv(file, sep=";")
            metadata["Sample number"] = metadata["Sample number"].astype(str)
            metadata = pd.merge(metadata, newsamples, on="Sample number", how="outer")
            metadata = metadata.rename(columns={"Spital: Ambulant/Stanionär": "treatment_type"})
            metadata.to_csv(args.outdir + "/" + plate + "/" + plate + "_metadata.csv", sep=";", index=None)

            not_empty_samples = [ find_ethid(sample[0], newsamples) for sample in all_samples_status if sample[1] == "not_empty"]
            #not_empty_anon = [ newsamples.loc[newsamples['Sample number'] == sample]["ethid"].to_string(index=False) for sample in not_empty_samples ]

            samplemap = create_sample_map(not_empty_samples, args.outdir)
            samplemap.rename(columns={0: "sample", 1: "lane"})
            try:
                samplemap.to_csv(args.outdir + "/" + plate + "/" + plate  + "_samplemap.tsv", sep="\t", index=None)
            except:
                sys.exit("Error: could not write the samplemap table. All directories and new lines in the pseudoanonymization table have been created and need to be deleted for the procedure to move forward.")
        else:
            print("Warning: all new samples appear to be already in the pseudo-anonymization table.")
            print("This is usually caused by samples that have been already pseudo-anonymized, but never processed by the pipeline.")
            print ("Skipping pseudo-anonymization")
    
