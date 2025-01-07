######
## Script name: gather_references.py
## Date: August 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Loads the list of
##        strains detected by Kraken2, checks the available references
##        and downloads missing ones
######

import pandas as pd
import argparse, os, sys, re, gzip, shutil, urllib
from Bio import SeqIO


# Script
if __name__ == '__main__':
	# Parse input args
    parser = argparse.ArgumentParser(description='Filters the full consensus list to keep only the consensus sequences of interest', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--substrains', required=True, type=str, help="The list of available substrains")
    parser.add_argument('--downloadlist', required=True, type=str, help="The list from the Kraken2 database containing the links for downloading the FASTA references")
    parser.add_argument('--fastadir', required=True, type=str, help="The directory for the fasta references")
    parser.add_argument('--bed', required=True, type=str, help="The directory for the bed references")
    parser.add_argument('--fastalinksdir', required=True, type=str, help="The sample directory for the links to the fasta references")
    args = parser.parse_args()
    
    #Some substrains are in the kraken2 database but not in the list of downloadable references from the same database
    #We skip these assignments
    stop_flag = 0
    with open("/data/config/missing.txt", 'r') as f:
        substrains_blacklist = f.read()
    substrains_blacklist = substrains_blacklist.split("\n")
    try:
        substrains = pd.read_table(args.substrains)
    except pd.errors.EmptyDataError:
        if os.path.exists(args.fastalinksdir) == False:
            os.mkdir(args.fastalinksdir)
        open(args.bed, 'a').close()
        print("WARNING: No substrains available. This is expected only in negative samples!")
        stop_flag = 1
    if stop_flag == 0:
        download_list = pd.read_table(args.downloadlist)        
        substrains = substrains.to_dict(orient="index")
        bed_list = []
        for i in range(0, len(substrains)):
            substrain = substrains[i]
            name = substrain['name,taxon_id'].split(",")[0]
            if name in substrains_blacklist:
                continue
            taxon_id = substrain['name,taxon_id'].split(",")[1]
            print("no fasta reference found for " + name + ". Searching if a download link is available")
            download_url = download_list[download_list["Sequence Name"].str.lower().str.contains(re.escape(name.lower()))]
            if len(download_url) > 1:
                test1 = download_url[download_url["Sequence Name"].str.lower().str.contains(re.escape(name.lower() + ","))]
                test2 = download_url[download_url["Sequence Name"].str.lower().str.contains(re.escape(name.lower() + " "))]
                if len(test1) == 1:
                    download_url = test1
                elif len(test2) ==1:
                    download_url = test2
                else:
                    continue
            download_url = download_url.reset_index()
            download_url = download_url.to_dict(orient="index")
            #for j in range(0, len(download_url)):
            try:
                this_url = download_url[0]
            except:
                print("ERROR!!!!" + name)
                with open("/data/config/missing.txt", 'a') as f:
                    f.write(name + "\n")
                sys.exit()
            print("Found the reference at URL: ", this_url["URL"])
            filename = os.path.join(args.fastalinksdir, str(taxon_id) + ".fasta.gz")
            filename_decompressed = re.sub('\.gz$', '', filename)
            if os.path.isfile(filename):
                os.remove(filename)
            urllib.request.urlretrieve(this_url["URL"], filename)
            with gzip.open(filename, 'rb') as f_in:
                with open(filename_decompressed, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            multifasta = SeqIO.parse(filename_decompressed, "fasta")
            #updated_multifasta = []
            for record in multifasta:
                bed_list.append({0:name, 1: 0, 2: len(record.seq), 3: record.description})
            #    record.description = ""
            #    updated_multifasta.append(record)
            #with open(filename_decompressed, "w") as output_handle:
            #    SeqIO.write(updated_multifasta, output_handle, "fasta")
            os.remove(os.path.join(args.fastalinksdir, filename_decompressed))
        bed = pd.DataFrame(bed_list)
        bed.to_csv(args.bed, index=None, header=None, sep="\t")

