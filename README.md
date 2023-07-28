# 2023_cb_stadler_ReVSeq
Code related to the project 2023_cb_stadler_ReVSeq.

## Introduction
This is a workflow design to perform sample retrieval and bioinformatics analysis for the project 2023_cb_stadler_ReVSeq. The goal of the project is to analyse sequencing data of patients' oropharingeal samples and detect which (if any) of the respiratory viruses of interest is the cause of infection observed by the physician.

The workflow is divided in two main sections:
- `pipeline_revseq`: the snakemake pipeline that performs the full bioinformatics analysis
- `revseq`: the docker container that retrieves samples and runs the pipeline

## Installation
### Pre-requisites
The workflow requires a working installation of docker. Please refer to the docker documentation for more information: https://docs.docker.com/engine/install/

### Installing and preparing the workflow

- Clone the repository
```
git clone git@github.com:ETH-NEXUS/2023_cb_stadler_ReVSeq.git
```
- Create the folders necessary to store the external data. It is suggested to create these folders outside of the git repository to avoid accidental commit of sensitive data or exposing secrets. There is no restriction for name and location because they are going to be added to the configuration file in a following step. Four folders are necessary:
  - `config`
  - `raw_data`
  - `metadata`
  - `anondata`
  - `viollier_mirror`
- Create the following four files and place them into the newly created folder `config`:
  - Create the anonymization table following the file `pipeline_revseq/config_templates/anonymization_table_template.tsv`. An example is provided in `pipeline_revseq/config_templates/anonymization_table_example.tsv`.
  - Create the sample map table following the file `pipeline_revseq/config_templates/sample_map_template.tsv`. An example is provided in `pipeline_revseq/config_templates/sample_map_example.tsv`
  - Create the snakemake config file following the file `pipeline_revseq/config_templates/config_example.yaml`.
  - Create an empty file called `empty_sample.tsv`.
- Create a file called `revseq/docker-compose.yml` following the file `revseq/docker-compose_example.yml` to match your system. The file must be in the git repository folder `revseq` and is already included in the `.gitignore` to avoid accidental commits 

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
- `uploadviollier`: Uploads back to Viollier the data the workflow is require to report back
- `anonymize_samples`: Anonymizes sample names by running the script `anonymize_sample_names.py`. Sample names are assigned a random, unique 6-characters ETHID from an alphabet composed of lowercase letters and numbers
- `dryrun`: Triggers a snakemake dry run of the bioinformatics pipeline. Useful for debugging and development
- `runpipeline`: Triggers a run of the bioinformatics pipeline
- `packagedata`: Packages data for import into the database

## Bioinformatics pipeline (pipeline_revseq folder)
A snakemake workflow is available to run the full analysis on the available raw data. It requires a `sample_map` file with a complete list of samples to analyse, an `anonymization_table` file with a list of all lanes associated to each sample, corretly configured `config.yaml`, `Dockerfile` and `docker-compose.yml` files.

The pipeline can be started at any moment using the command `./revseq/revseq runpipeline`.

### Available steps
The pipeline is divided in the four major areas: `preprocessing`, `basic pipeline`, `dehumanization`, `qc`, with the following steps:
-  `preprocessing`: prepare the raw data for the bioinformatics analysis
  - `merge_lanes`: merge the raw data from multiple lanes in a single fastq file
  - `bwa_index`: index the virus reference sequences for use with BWA
  - `trim_galore`: trim primers from the raw sequences and filter reads by quality and by read length
- `basic pipeline`: process the samples to detect the virus(es) present
  - `bwa`: align the trimmed reads to the virus refernce sequences using BWA
  - `remove_multimappers`: remove reads mapping in multiple contigs (i.e. multiple viruses from the alignment results)
  - `sort`: Sort the alignments by position
  - `remove_duplicate`: remove PCR and optical duplicates from the alignment results
  - `samtools_index`: create the BAM index for the deduplicated alignments
  - `pileup`: compute coverage information on viral sequences with samtools mpileup for the deduplicated alignments
  - `idxstats`: compute and report index stats and aligned reads counts for the deduplicated alignments
  - `assign_virus`: assign the virus(es) more likely to be present in the sample. More details on the custom method available in the next paragraph
  - `validate_assignment`: compare the virus assignment with the results available from the screening panels and add the validation to the report
- `dehumanization`: remove human reads from the raw data for publication
  - `dh_reuse_alignreject`: extract the unaligned reads from the alignment results
  - `dh_host_index`: index the host reference genome for use with BWA
  - `dh_hostalign`: align to the host genome the reads that did not align to the virus sequences
  - `dh_filter`: filter from the raw FASTQ files the reads that align on the host genome
  - `dehuman`: package the filtered reads in cram format
- `qc`: Quality control steps and report
  - `fastqc`: qaulity control of raw reads using FastQC
  - `samtoolsstats`: quality control of the deduplicated alignments using samtoolsstats
  - `rseqc`: quality control of the deduplicated alignments using rseqc
  - `qualimap`: quality control of the deduplicated alignments using qualimap
  - `multiqc`: collecting and reporting all previous QC steps using multiQC
  - `fastqc_dehuman`: quality control of the dehumanized FASTQ files using FastQC
  - 
A rulegraph showing how these steps are connected is available in `imagesrevseq_pipeline_dag.pdf`

#### Details - Virus assignment
The virus assignment is computed by a custom python script that summarises the results in a dedicated table.
From the index stats, the tool retrieves the amount of reads aligned to each virus sequence and the length of each viral sequence.
Virus sequences related to different ORF of the same virus strain are merged together by sum of the all the ORF's lengths and sequences.
This information is then used to compute RPKMs and the percentage of RPKM assigned to each virus sequence.

Virus sequences with outlier RPKM values are detected by searching for RPKM percentage larger than a user-defined higher percentile of the RPKM percentage distribution. The default is the 90th percentile.

The results are saved in a dedicated table, where outliers are marked with a `*`.

If the user requested it, strains are then collapsed together (e.g. all strains referring to different sub-strains of Influenza A can be collapsed to a single cumulative values for Influenza A in general). RPKM, RPKM percentages and outliers are re-computed based on the new read counts and lengths.

The collapsed results are then compared with the results available in the respiratory virus test panels run on the samples by the clinicians. Strains reported as positive in the test panels are reported with a `*`.

# Acknowledgements
- [SVVC - Neher Lab](https://github.com/neherlab/SVVC): for providing the base concept on which this method has been developed
- [V-Pipe - ETHZ CBG](https://github.com/cbg-ethz/V-pipe): for providing the rules for dehumanization
