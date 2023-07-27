######
## Script name: link_restructure.py
## Date: July 2023
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq command that controls the revseq workflow.
##      Takes the Viollier mirror samples and restuctures it using softlinks
##      in a new directory to match the pipeline requirements
######

import sys, os
import argparse

# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference',
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--mirror', required=True, type=str, help='the directory containing the Viollier mirror')
    parser.add_argument('--sampledir', required=True, type=str, help='the directory containing all raw data samples, one folder per sample')

    args = parser.parse_args()
    
    raw_data = [ file.name for file in os.scandir(args.mirror) ]
    linked_files = []
    linked_samplenames = []
    for folder in os.scandir(args.sampledir):
        linked_samplenames=linked_samplenames.append(folder.name)
        for files in os.scandir(args.sampledir + "/" + folder.name):
            linked_files=linked_files.append(files.name)

    to_link = list(set(raw_data).difference(linked_files))

    for file in to_link:
        target_dir = file.split("_")[0]
        if not os.path.exists(sampledir + "/" + target_dir):
            os.mkdir(sampledir + "/" + target_dir)
        else:
            if not os.path.isdir(sampledir + "/" + target_dir):
                sys.exit("ERROR: tried to create directory " + target_dir + " but a file with the same name already exists")
                break
        if os.path.exists(sampledir + "/" + target_dir + "/" + file):
            sys.exit("ERROR: tried to link an already existing file: " + file + ". This should NEVER happen as this is checked beforehand!")
        try:
            os.symlink(mirror + "/" + file, sampledir + "/" + target_dir + "/" + file)
        except:
            sys.exit("ERROR: filed to link sample " + file)
        print("Linked file " + file)

