# 2023_cb_stadler_ReVSeq
Code related to the project 2023_cb_stadler_ReVSeq.

## Introduction
This is a workflow design to perform sample retrieval and bioinformatics analysis for the project 2023_cb_stadler_ReVSeq. The goal of the project is to analyse sequencing data of patients' oropharingeal samples and detect which (if any) of the respiratory viruses of interest is the cause of infection observed by the physician.

The workflow is divided in two main sections:
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
- The workflow works with data that may be sensitive, e.g. paths, internal sample names, metadata. It is good practice to leave this data outside of the git repository to avoid accidental commits. We provide a script to simplify the procedure
  - Run the script `resoruces/initialize.sh` and follow the instructions
    - The script creates and populates in the provided directory the following 5 sub-directories:
      - `pipeline_configuration`: where the pipeline configuration is stored
        - The file created in `pipeline_configuration` are already customized based on user input
      - `viollier_mirror`: where raw data and metadata from the Viollier SFTP server will be mirrored
      - `raw_data`: where the mirrored data will be re-organized to be more user- and pipeline-friendly
      - `anondata`: where the samples will be anonymized. This will create a physical separation between non-anonymized and anonymized data
      - `results`: where the bioinformatics pipeline analysis results will be written
    - The script also creates and customizes the content of the file `revseq/docker-compose.yml` based on the user input

### Running the workflow
Run the container by navigating to folder `revseq` and running the command
```
docker-compose up --detach
```

## Revseq command (revseq folder)
The `revseq` command is the hub that controls all workflow steps. When the container is run, the main commands are called automatically by the script `revseq/run_revseq.sh`, but the `revseq` command can also be used as a standalone tool to call specific workflow step.

### Available sub-commands
- `syncviollier`: calls lftp to connect to Viollier's sftp server and mirror the directories where the data for the Revseq project are stored.
- `links_amples`: The pipeline requires a specific raw data and metadata directory structure to function properly. This command uses the script `link_restructure.py` to take the samples in the Viollier mirror and create the necessary directory structure in the `links` folder by using softlinks
- `uploadviollier`: Uploads data  back to Viollier that the workflow is required to report back
- `anonymize_samples`: Anonymizes sample names by running the script `anonymize_sample_names.py`. Sample names are assigned a random, unique 6-characters ETHID from an alphabet composed of lowercase letters and numbers
- `dryrun`: Triggers a snakemake dry run of the bioinformatics pipeline. Useful for debugging and development
- `runpipeline`: Triggers a run of the bioinformatics pipeline

## Bioinformatics pipeline (pipeline_revseq folder)
A snakemake workflow is available to run the full analysis on the available raw data. It requires a `sample_map` file with a complete list of samples to analyse, an `anonymization_table` file with a list of all lanes associated to each sample, corretly configured `config.yaml`, `Dockerfile` and `docker-compose.yml` files.

The pipeline can be started at any moment using the command `./revseq/revseq runpipeline`.

### Available steps
The pipeline is divided in the four major areas: `preprocessing`, `basic pipeline`, `dehumanization`, `qc`, with the following steps:
-  `preprocessing`: prepare the raw data for the bioinformatics analysis
  - `merge_lanes`: merge the raw data from multiple lanes in a single fastq file
  - `bwa_index`: index the merged host-virus reference sequences for use with BWA
  - `trim_galore`: trim primers from the raw sequences and filter reads by quality (default) and by read length (minimum 80bp after trimming)
- `dehumanization`: remove human reads from the dataset
  - `bwa`: align the trimmed reads to the merged host-virus referenece sequences using BWA
  - `dh_postprocess`: sort and index alignments
  - `dehuman`: remove any read wit primary, secondary or multiple alignments on any host sequence
  - `bam_to_fastq`: convert the dehumanized BAM into FASTQ for CRAM compression
  - `cram`: compressed the dehumanized FASTQs in CRAM format
- `basic pipeline`: process the samples to detect the virus(es) present
  - `filter_alignment`: remove reads with secondary, multiple or split alignments. Remove reads with quality 0
  - `sort`: Sort the alignments by position
  - `remove_duplicates`: remove PCR and optical duplicates from the alignment results
  - `samtools_index`: create the BAM index for the deduplicated alignments
  - `pileup`: compute coverage information on viral sequences with samtools mpileup for the deduplicated alignments
  - `idxstats`: compute and report index stats and aligned reads counts for the deduplicated alignments
  - `assign_virus`: assign the virus(es) more likely to be present in the sample. More details on the custom method available in the next paragraph
  - `validate_assignment`: compare the virus assignment with the results available from the screening panels and add the validation to the report
- `qc`: Quality control steps and report
  - `fastq_raw`: quality control of raw reads using FastQC
  - `fastqc_merged`: quality control of lane-merged reads using FastQC
  - `dh_fastqc`: quality control of dehumanized reads using FastQC
  - `samtoolsstats`: quality control of the deduplicated alignments using samtoolsstats
  - `rseqc`: quality control of the deduplicated alignments using rseqc
  - `qualimap`: quality control of the deduplicated alignments using qualimap
  - `multiqc`: collecting and reporting all previous QC steps using multiQC
  
A rulegraph showing how these steps are connected is available in `images/revseq_pipeline_rulegraph.pdf`

#### Details - Virus assignment
The virus assignment is computed by a custom python script that summarises the results in a dedicated table.
From the index stats, the tool retrieves the amount of reads aligned to each virus sequence and the length of each viral sequence.
Virus sequences related to different ORF of the same virus strain are merged together by sum of the all the ORF's lengths and sequences.
This information is then used to compute RPKMs and the percentage of RPKM assigned to each virus sequence.

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