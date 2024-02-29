######
## Script name: upload_to_db.py
## Date: February 2024
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Part of the Revseq pipeline. Uses the functions available in
##        RevSeqDataUploader (https://github.com/ETH-NEXUS/RevSeqDataLoader)
##        to upload a plate's results to the database
######

import argparse, os, sys

# Script
if __name__ == '__main__':
    # Parse input args
    parser = argparse.ArgumentParser(description='fetch the primers positions on the reference', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--plate', required=True, type=str, help='Barcode of the plate to upload')
    parser.add_argument('--revseqdataloader', type=str, default='/app/RevSeqDataLoader', help='Path of the directory containing the RevSeqDataLoader package')
    args = parser.parse_args()


sys.path.append(args.revseqdataloader)
from RevSeqDataLoader import RevSeqDataLoader

# Initialize the RevSeqDataLoader
loader = RevSeqDataLoader()

# Specify the path to the data you want to import
path = "/data/" + args.plate

# Call the import_data method with the path
try:
    loader.import_data(path)
except:
    sys.exit("Error in loading data to the database")
    