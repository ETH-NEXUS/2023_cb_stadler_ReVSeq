# 2023_cb_stadler_ReVSeq
Code related to the project 2023_cb_stadler_ReVSeq

# Snakemake workflow
A snakemake workflow is available to run the full analysis on the available raw data. It requires a `sample_map` file with a complete list of samples to analyse and can be started at any moment using the command `./revseq/revseq runpipeline`.

## Available steps
The pipeline is divided in the following steps:
- `fetch_primers`: retrieve the positions of the primers on the reference genome
- `bwa_index`: index the reference genome
- `trim_galore`: quality and length trimming of raw reads
- `bwa`: align reads on reference genome
- `fasta_index`: use samtools to index the primers' fasta file
- `pileup`: use samtools to retrieve allele counts
- `pair_frequences`: use a custom script to retrieve the frequencies
- `consensus`: use a custom script to create a consensus sequence
- `dehuman`: filter out human reads from the raw data for publication

# Acknowledgements
- [SVVC - Neher Lab](https://github.com/neherlab/SVVC): for providing the base on which this method has been developed
- [V-Pipe - ETHZ CBG](https://github.com/cbg-ethz/V-pipe): for providing the rules for filtering out human reads