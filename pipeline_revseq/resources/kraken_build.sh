#!/bin/bash
######
## Script name: kraken_build.sh
## Author: Matteo Carrara (NEXUS Personalized Health Technologies)
##
## Description: Builds a custom Kraken2 database starting from a list of 
##        FASTA references
######

# Collect files into an array (modify the glob pattern as needed)
conda activate kraken2
files=(./*fasta)
total=${#files[@]}
counter=0

# Variable to track the last percentage that was printed
last_printed=0

for file in "${files[@]}"; do
    # Place your file-processing logic here
    # e.g., process_file "$file"
    kraken2-build --add-to-library "$file" --db revseq_custom
    ((counter++))
    percent=$(( counter * 100 / total ))

    # If the percent is a new multiple of 5, print the progress
    if (( percent >= last_printed + 5 )); then
        last_printed=$(( (percent / 5) * 5 ))
        printf "\rProcessing file %d of %d (%d%%)" "$counter" "$total" "$last_printed"
    fi
done

# Ensure we print 100% if not already printed, and then move to a new line
if (( last_printed < 100 )); then
    printf "\rProcessing file %d of %d (100%%)" "$total" "$total"
fi
echo -e "\nDone!"

echo "Building the database"
kraken2-build --download-taxonomy --db revseq_custom
kraken2-build --build --db revseq_custom

