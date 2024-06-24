# Content:
The file cds.bed contains the CDS annotations for all viruses included in the ReVSeq project.

# Base annotation retrieval
Data has been downloaded from NCBI using the command line tool "datasets" version 14.3.2.
The date of download is June 5th 2024.
Using all available tax IDs, the data has been retrieved with the command

```
datasets download virus genome taxon <taxon> --include annotation --filename <filename>.zip
```
After de-compression, only the files named `data/annotation_report.jsonl` have been kept.

### Exceptions
- It was not possible to download record NC_004148.2 because it was removed from RefSeq.
- It was not possible to extract the CDS information from the records NC_011203, NC_009238, CY115156.1, NC_003266, NC_001405, NC_001617, NC_001803.
- Record MN908947.3 was too large
In those specific cases, the CDS annotation has been manually retrieved from the NCBI website

# Downstream
From the jsonl files, we are interested only in the CDS positions. Script find_cds.py retrieves
the CDS positions for any jsonl file.
