# ReVSeq - Respiratory Virus Sequencing

## Introduction
This repository contains code for data retrieval and bioinformatics analysis of samples from patients with respiratory virus symptoms tested and sequenced at Viollier AG. The project aims at providing an overview on the genetic diversity of respiratory viruses in Switzerland as well as fostering an understanding of their evolution.

The analysis workflow is divided in two main sections:
- `pipeline_revseq`: the snakemake pipeline that performs the full bioinformatics analysis
- `revseq`: the docker container that retrieves samples and runs the pipeline

## Installation
### Pre-requisites
- The workflow requires a working installation of docker. Please refer to the docker documentation for more information: https://docs.docker.com/engine/install/
- The worklow requires the configuration necessary to connect to the Viollier SFTP server. This includes
  - `.ssh/config` file with a section dedicated to the connection to the SFTP server
  - Private and public key pair to be used for connection. For the pair to be usable, the public key must be shared with Viollier
  - `.ssh/known_hosts` file with a section dedicated to the Viollier SFTP server public key

### Installing and preparing the workflow

- Clone the repository
```
git clone git@github.com:ETH-NEXUS/2023_cb_stadler_ReVSeq.git
```
- Run the script `resources/initialize.sh` and follow the instructions
    - The script creates and populates in the provided directory the following 5 sub-directories:
      - `pipeline_configuration`: where the pipeline configuration is stored
        - The file created in `pipeline_configuration` are already customized based on user input
      - `viollier_mirror`: where raw data and metadata from the Viollier SFTP server will be mirrored
      - `raw_data`: where the mirrored data will be re-organized to be more user- and pipeline-friendly
      - `anondata`: where the samples will be anonymized. This will create a physical separation between non-anonymized and anonymized data
      - `results`: where the bioinformatics pipeline analysis results will be written
    - The script also creates and customizes the content of the file `revseq/docker-compose.yml` based on the user input

### Running the full workflow
By default, the workflow is set to run within a docker container, where a monitoring loop regularly updates the mirror of the raw data from the sequencing center, pseduonymizes sample names, and runs the full bioinformatics pipeline.

To run the bioinformatics pipeline independently from raw data mirroring, please refer to the section `Running the bioinformatics pipeline manually` below.


Run the container by navigating to the folder `<git repo>/revseq` and running the command
```
docker-compose up --detach
```
to start the docker container (default name `revseq`) responsible to run the workflow.

#### Revseq container
The revseq docker container starts by running an infinite loop responsible of monitoring the raw sequencing data mirror, pseduonymize samples, and run the snakemake pipeline on the any newly detected sample.

##### Commands
The loop takes care of running automatically the following three commands:
- `<git repo>/revseq/revseq syncvioller`: checks on the sequencing center server if new raw sequencing data is available and mirrors it to the local folder defined in the configuration
- `<git repo>/revseq/revseq pseudoanonymize_samples`: assigned a pseudonymized ID to each new sample. The pseudonym is ensured to be unique and is an alphanumeric string of 6 characters that can include both uppercase and lowercase letters
- `<git repo>/revseq/revseq runpipeline \<plate\>`: runs the full bioinformatics pipeline on the all pseudonymized samples included in the delivery batch

### Running the bioinformatics pipeline manually
The bioinformatics analysis step is not strictly dependent on a direct connection to the seqencing center, nor on the pseudonymization of the samples. The command `<git repo>/revseq/revseq revseq \<plate\>` can in fact be run manually to trigger only the bioinformatics analysis outside of the monitoring loop. This section explains how to start a docker container that can be used for such purpose.

Run the container by navigating to the folder `<git repo>/revseq/manual_run` and running the command
```
docker-compose up --detach
```
to start the docker container (default name `revseq_manual`) responsible to wait for user input and configured to run the bioinformatics pipeline.

#### revseq_manual container
The revseq_manual docker container starts idle, configured to run any command from `<git repo>/revseq/manual_run/revseq (`/data/revseq` inside the container).

##### How to run the commands
With the docker container `revseq_manual` up and running, open an interactive session within it by running the command
```
docker exec -it revseq-revseq_manual-1 bash
```
The new prompt is within the container and allows access to its full configuration.

To run the pipeline it is only necessary to run:
```
/data/revseq/revseq pseudoanonymize_samples
/data/revseq/revseq runpipeline <batch>
```

## Input
The workflow requires several input files:
- raw FASTQ data: raw sequencing reads in FASTQ format. The pipeline only accepts paired-end protocols and expects 2 files per samples
- \<batch\>\_samplemap.tsv: a tab-separated table with two columns and no header. The first column lists the sample names and the second lists one of the available lane codes for the samples. The file guides the merging of multi-lane raw data


## Output
The workflow analyses the raw sequences data and provides several output files that are part of the quantification and characterization the viruses detected in each sample. 
The workflow provides outputs for all analysis steps. The final outputs of interest are:
- \<batch\>\_empty_samples.txt: a list of all samples found with zero raw reads. Such samples are not analyzed by the bioinformatics workflow even if they are included in the sample list.
- \<batch\>\_pseudoanon_table.tsv: a table matching the original sample ID to its pseudonymized ID
- \<batch\>\_metadata.csv: a copy of the clinical metadata including the pseudonymized ID
- gather_results: directory gathering all main results in one place
  - sample\_\<samplename\>: directory containing the results for the specific sample _samplename_
    - chr_file.txt: file listing all available chromosome of the major assigned virus
    - \<sample\>\_consensus.fa: fasta file with the consensus sequence of all major assigned viruses
    - \<sample\>\_consensus_upload.gz: compressed fasta file with the consensus sequence of the major assigned virus
    - \<sample\>\_count_n.txt: table listing the number and proportion of Ns (masked positions due to low coverage) in the consensus sequences of all major assigned viruses
    - \<sample\>\_count_table.tsv: full virus assignment table listing the assignment results for all panel viruses plus the viruses detected by Kraken2
    - \<sample\>.cram: cram file containing the dehumanized raw reads
    - \<sample\>\_merged_R[1|2].fastq.gz: compressed FASTQ files containing the lane-merged forward and reverse raw reads
    - \<sample\>\_remove_duplicated.bam: final bam file containing the aligned reads after cleanup and filtering


## Revseq command (revseq folder)
The `revseq` command is the hub that controls all workflow steps. When the container is run, the main commands are called automatically by the script `revseq/run_revseq.sh`, but the `revseq` command can also be used as a standalone tool to call specific workflow steps.

### Available sub-commands
- `syncviollier`: calls lftp to connect to Viollier's sftp server and mirror the directories where the data for the Revseq project are stored.
- `link_samples`: The pipeline requires a specific raw data and metadata directory structure to function properly. This command uses the script `link_restructure.py` to take the samples in the Viollier mirror and create the necessary directory structure in the `links` folder by using softlinks
- `uploadviollier`: Uploads required results to Viollier
- `anonymize_samples`: Anonymizes sample names by running the script `anonymize_sample_names.py`. Sample names are assigned a random, unique 6-characters ETHID from an alphabet composed of lowercase letters and numbers
- `dryrun`: Triggers a snakemake dry run of the bioinformatics pipeline. Useful for debugging and development
- `runpipeline`: Triggers a run of the bioinformatics pipeline
- `unlock`: Manually unlocks the working directory
- `list_temp_files`: Lists the temporary files generated by the pipeline
- `delete_temp_files`: Deletes the temporary files generated by the pipeline

## Bioinformatics pipeline (pipeline_revseq folder)
A snakemake workflow is available to run the full analysis on the available raw data. It requires a `sample_map` file with a complete list of samples to analyse, an `anonymization_table` file with a list of all lanes associated to each sample, correctly configured `config.yaml`, `Dockerfile` and `docker-compose.yml` files.

The pipeline can be started at any moment using the command `./revseq/revseq runpipeline`.

### Available steps
![RuleGraph](images/revseq_pipeline_rulegraph.png?raw=true "Rulegraph")
The pipeline is divided in the five major areas: `preprocessing`, `dehumanization`, `basic_pipeline`, `qc`, `aggregation`, and `uploads` with the following steps:
-  `preprocessing`: prepare the raw data for the bioinformatics analysis
  - `merge_lanes`: merge the raw data from multiple lanes in a single fastq file
  - `kraken2`: metagenomic strain assignment. The reported strains are used downstream to prepare the reference for the alignment
  - `gather_references`: generate a bed file of all the viruses detected by kraken2
  - `merge_refs`: generate a multi-fasta reference file of all the viruses detected by kraken2
  - `bwa_index`: index the merged host-virus reference sequences for use with BWA
  - `cutadapt`: trim primers and adapters from the raw sequences and filter reads by quality (default) and by read length (minimum 80bp after trimming). Adapters are trimmed based on an optional metadata file. If the file is absent, no adapter trimming is performed
- `dehumanization`: remove human reads from the dataset
  - `bwa`: align the trimmed reads to the merged host-virus referenece sequences using BWA
  - `dehuman`: remove any read wit primary, secondary or multiple alignments on any host sequence
  - `bam_to_fastq`: convert the dehumanized BAM into FASTQ for CRAM compression
  - `cram`: compressed the dehumanized FASTQs in CRAM format
- `basic pipeline`: process the samples to detect the virus(es) present
  - `filter_alignment`: remove reads with secondary, multiple or split alignments. Remove reads with quality 0
  - `remove_duplicates`: remove PCR and optical duplicates from the alignment results
  - `merge_regions`: postprocess the output of gather_references to follow the expected format
  - `samtools_index`: create the BAM index for the deduplicated alignments
  - `pileup`: compute coverage information on viral sequences with samtools mpileup for the deduplicated alignments
  - `depth`: compute the depth information on viral sequences with samtools depth for the deduplicated alignments
  - `idxstats`: compute and report index stats and aligned reads counts for the deduplicated alignments
  - `assign_virus`: assign the virus(es) more likely to be present in the sample. More details on the custom method available in the next paragraph
  - `validate_assignment`: compare the virus assignment with the results available from the screening panels and add the validation to the report
  - `consensus`: create the consensus sequence for all available references using BCFTools. Any position with coverage below the minimum coverage defined in the configuration yaml file is masked with a 'N' character
  - `postprocess_consensus`: filter the consensus to keep only the references of interest. Count the statistics of the consensus (number of Ns, proportion of Ns, coverage of non N regions)
  - `consensus_cds`: create the consensus sequence for all available references of only the CDS regions. The consensus is generated starting from the results of the step `consensus` and filtering based on the BED file listing all CDS regions. The result is a multiFASTA listing one CDS per entry
  - `count_n_cds`: Count the statistics of the CDS-based consensus (number of Ns, proportion of Ns, coverage of non N regions)
  - `chr_file`: Generate the chromosome file necessary for the ENA uploads
- `qc`: Quality control steps and report
  - `fastq_raw`: quality control of raw reads using FastQC
  - `fastqc_merged`: quality control of lane-merged reads using FastQC
  - `fastqc_trimmed`: quality control of primer-trimmed reads using FastQC
  - `dh_fastqc`: quality control of dehumanized reads using FastQC
  - `samtoolsstats`: quality control of the deduplicated alignments using samtoolsstats
  - `samtoolsstats_filtered`: quality control of the filtered alignments using samtoolsstats
  - `rseqc`: quality control of the deduplicated alignments using rseqc
  - `rseqc_filtered`: quality control of the filtered alignments using rseqc
  - `qualimap`: quality control of the deduplicated alignments using qualimap
  - `qualimap_filtered`: quality control of the filtered alignments using qualimap
  - `multiqc`: collecting and reporting all previous QC steps using multiQC
  - `multiqc_filtered`: collecting and reporting all previous QC steps for the filtered alignments using multiQC
- `aggregation`: rules to aggregate the plate results for a human-readable summary of assignments
  - `aggregate`: generate the tables summarizing the virus assignment and the qc
- `upload`: functions to upload results on the database and on Viollier's FTP server
  - `gather_results_samples`: processes the assignment table to add the consensus statistics. Gather all sample-related files for upload in a single directory
  - `gather_results_plate`: Gather all plate-related files for upload in a single directory
  - `push_to_db`: use ReVSeqDataLoader to upload the gathered data to the database
  
A rulegraph showing how these steps are connected is available in `images/revseq_pipeline_rulegraph.pdf`

#### Details - Kraken2
The Kraken2 assignment is performed using a curated custom database generated as follows:
- Downloaded the Kraken2 standard reference https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz
- From the full set of references we removed:
  - Archaea, Plasmids, all animals except for human (reason: not part of the project)
  - STI-related viruses from https://www.niaid.nih.gov/diseases-conditions/sti-pathogens-and-syndromes and Polyovirus (reason: potentially sensitive information)
  - Pandoraviruses (reason: not part of the project. The large genome increases the memory requirements unnecessarily)
  - Tupanvirus deep ocean, Prymnesium kappa virus, Cotonvirus japonicus, Bodo saltans virus, Tupanvirus soda lake, Powai lake megavirus, Acanthamoeba polyphaga mimivirus (reason: broken download links)
  - Human coronavirus OC43, Severe acute respiratory syndrome coronavirus 2, Respiratory syncytial virus (reason: already in the Illumina panel and always included in the analysis)
  - Simian Agent 10, Bovine coronavirus (reason: too similar in sequence to viruses in the Illumina panel)
- We use the metadata file `resources/library_report.tsv` to download the FASTA reference of all the species we kept in our custom list
- We use the custom script `resources/kraken_build.sh` to add all reference FASTA sequences to a new database and then build it using kraken-build

#### Details - Virus assignment
The virus assignment is computed by a the python script `pipeline_revseq/workflow/scripts/assign_virus.py` that summarises the results in a dedicated table.
From the index stats, the tool retrieves the amount of reads aligned to each virus sequence and the length of each viral sequence.
Virus sequences related to different Open Reading Frames (ORF) of the same virus strain are merged together by sum of the all the ORF's lengths and sequences.
This information is then used to compute Reads per Kilobase per Million mapped reads (RPKMs) and the percentage of RPKM assigned to each virus sequence.

Virus sequences with outlier RPKM values are detected by searching for RPKM percentage larger than a user-defined higher percentile of the RPKM percentage distribution. The default is the 90th percentile.

The results are saved in a dedicated table, where outliers are marked with a `*`.

If the user requested it, strains are then collapsed together (e.g. all strains referring to different sub-strains of Influenza A can be collapsed to a single cumulative values for Influenza A in general). RPKM, RPKM percentages and outliers are re-computed based on the new read counts and lengths.

The collapsed results are then compared with the results available in the respiratory virus test panels run on the samples by the clinicians. Strains reported as positive in the test panels are reported with a `*`.

# Web App Revseq: Integration of API and UI for Database Output File Imports 

## Initial Deployment 

Upon the initial deployment of the application, certain static content must be populated into the database. This includes the data for:

- Panel-Strain Lookup 
- Strain-Substrain Lookup 
- Strain-ORF Lookup 
- File Type Postfixes 

To load this static content into the database, execute the commands below within the API Docker container:

```bash
python manage.py fixture panel_strain
python manage.py fixture strain_substrain
python manage.py fixture file_type
```

## Importing Data into the Database 

Data import is organized by plates. The command for importing requires a single argument: the directory path containing the results, which is typically named after the plates' barcode. 

```bash
python manage.py import <path_to_results_dir>
```


## Authentication & Login

### Browser Access

When accessing the API through a browser, you'll automatically be redirected to the login page if you're not authenticated.

### Terminal Access

For terminal or command-line access, you'll need to use JWT (JSON Web Token) authentication. Follow the steps below:

1 Use the following command to obtain an access token:

```bash
curl -X POST -H "Content-Type: application/json"  -H "Accept: application/json" -d '{"username": "your_username", "password": "your_password"}' https://revseq.nexus.ethz.ch/api/token/

```

2. Include the obtained access token in the Authorization header of your requests. For example:

```bash
curl -X GET  -H "Authorization: Bearer your_access_token" -H "Accept: application/json"  https://revseq.nexus.ethz.ch/api/samples/


```

## Uploading Data to the Database

```
curl -X POST  -H "Authorization: Bearer your_access_token"  -H "Accept: application/json"   -H "Content-Type: application/json" -d '{"path":"/path/to/results/dir"}' https://revseq.nexus.ethz.ch/api/import_results/

```

## Change mode

This script is meant to be executed once. Enter the api container and run the following command:

```bash
python manage.py change_mode <mode> --add_prefix --prefix <prefix>
```
<mode>: The mode to set for samples. Must be either metagenomics or alignment.
--add_prefix: Optional flag to indicate whether to add a prefix to the pseudonymized_id.
--prefix <prefix>: The prefix to add to the pseudonymized_id if --add_prefix is set.

 docker compose  -f docker-compose.yml -f  docker-compose.dev.yml  up
