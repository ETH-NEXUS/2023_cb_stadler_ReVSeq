# Content:
The file cds.bed contains the CDS annotations for all viruses included in the ReVSeq project.

# Base annotation retrieval
Data has been downloaded from NCBI by navigating to the `nuccore` page for the specific virus ID (e.g. `https://www.ncbi.nlm.nih.gov/nuccore/NC_021928`) and using the download function (i.e. the section `send to` in the top right corner) to download the `complete record` in GFF3 format.
The date of download is June 24th 2024.

# Downstream
From the GFF3 files, we are interested only in the CDS positions to create a flat bed file.
First, BEDOPS 2.4.41 has been used to generate 0-based bed files from the GFF3. Then, grep is used to keep only lines describing a CDS. Bedtools 2.31.1 is used to merge the overlapping CDS entries to flatten the bed files. Finally, all bed files are merged in one reference bed.
The only exceptions are NC_001490, NC_001617, NC_009996, AB686524, GQ865517, for which we also need to keep the position of the 5'UTR regions, because the 5'UTR region is often used for identification.
Those strains are Rhino- and Entrero viruses. However, only the Enteroviruses (AB686524, GQ865517) included annotation for the 5'UTR.
