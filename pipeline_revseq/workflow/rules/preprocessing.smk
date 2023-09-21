## Helper functions
import os
import pandas as pd

def find_input_fastq_lanes_r1(wildcards):
    #"return a list of all input fastq files available for one sample. Handles multi-lane data"
    mate = "R1"
    sample = wildcards.sample
    matches = pd.read_table(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv")
    snumber = int(matches[matches["ethid"] == sample]["Sample number"].to_string(index=False))
    sample_table = sample_map[sample_map["sample"] == sample]
    sample_lanes = sample_table["lane"].tolist()
    sample_lanes.sort()
    fastq_files = []
    total_samples = os.listdir(config["inputOutput"]["input_fastqs"])
    for lane in sample_lanes:
        all_samples = [ file for file in total_samples if str(snumber) in file ]
        lanefile = [ file for file in all_samples if lane in file ]
        raw = [ file for file in lanefile if mate in file ]
        if len(raw) == 0:
            sys.exit("ERROR: can't file raw file for sample " + snumber + " lane " + lane + " mate " + mate)
        fastq_files.append(config["inputOutput"]["input_fastqs"]+"/"+raw[0])
    return fastq_files


def find_input_fastq_lanes_r2(wildcards):
    #"return a list of all input fastq files available for one sample. Handles multi-lane data"
    mate = "R2"
    sample = wildcards.sample
    matches = pd.read_table(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv")
    snumber = int(matches[matches["ethid"] == sample]["Sample number"].to_string(index=False))
    sample_table = sample_map[sample_map["sample"] == sample]
    sample_lanes = sample_table["lane"].tolist()
    sample_lanes.sort()
    fastq_files = []
    total_samples = os.listdir(config["inputOutput"]["input_fastqs"])
    for lane in sample_lanes:
        all_samples = [ file for file in total_samples if str(snumber) in file ]
        lanefile = [ file for file in all_samples if lane in file ]
        raw = [ file for file in lanefile if mate in file ]
        if len(raw) == 0:
            sys.exit("ERROR: can't file raw file for sample " + snumber + " lane " + lane + " mate " + mate)
        fastq_files.append(config["inputOutput"]["input_fastqs"]+"/"+raw[0])
    return fastq_files


## Rules
rule merge_lanes:
    input:
        r1 =find_input_fastq_lanes_r1,
        r2 =find_input_fastq_lanes_r2,
    output:
        r1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_lanes/{sample}_merged_R1.fastq.gz"),
        r2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_lanes/{sample}_merged_R2.fastq.gz"),
    params:
        rawdir = config["inputOutput"]["input_fastqs"],
        pseudoanon_table = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv",
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_lanes",
        outsuffix_r1 = "_merged_R1.fastq.gz",
        outsuffix_r2 = "_merged_R2.fastq.gz"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/merge_lanes/merge_lanes.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/merge_lanes/merge_lanes.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/merge_lanes/{sample}_merge_lanes.benchmark"
    shell:
        """
        python workflow/scripts/merge_lanes.py \
        --rawdir {params.rawdir} \
        --lanefiles {input.r1} \
        --pseudoanon_table {params.pseudoanon_table} \
        --outdir {params.outdir} \
        --outsuffix {params.outsuffix_r1}

        python workflow/scripts/merge_lanes.py \
        --rawdir {params.rawdir} \
        --lanefiles {input.r2} \
        --pseudoanon_table {params.pseudoanon_table} \
        --outdir {params.outdir} \
        --outsuffix {params.outsuffix_r2}
        """

rule merge_refs:
    input:
        virus_reference = config["resources"]["reference"],
        host_reference = config["resources"]["host_ref"],
    output:
        referenceout = config["resources"]["reference_dir"]+"/merged_virus_host_ref.fa",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/merge_refs/merge_refs.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/merge_refs/merge_refs.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/merge_refs/merge_refs.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        cat {input.virus_reference} > {output.referenceout} 2> >(tee {log.errfile} >&2)
        cat {input.host_reference} >> {output.referenceout} 2> >(tee {log.errfile} >&2)
        """


rule bwa_index:
    input:
        reference = config["resources"]["host_ref"]
    output:
        ref_index = config["resources"]["host_ref"] + '.bwt'
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/bwa_index/bwa_index.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/bwa_index/bwa_index.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/bwa_index/bwa_index.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.reference} 2> >(tee {log.errfile} >&2)"


rule trim_galore:
    input:
        r1 = rules.merge_lanes.output.r1,
        r2 = rules.merge_lanes.output.r2,
    output:
        r1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_galore/{sample}_merged_R1_val_1.fq.gz"),#temp
        r2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_galore/{sample}_merged_R2_val_2.fq.gz"),#temp
    params:
        base_name = "data/{sample}",
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_galore",
        min_length = config["tools"]["trim_galore"]["min_length"],
        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/trim_galore/{sample}_trim_galore.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/trim_galore/{sample}_trim_galore.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/trim_galore/{sample}.benchmark"
    conda:
        "../envs/trim_galore.yaml"
    threads: config['threads']
    shell:
        """
        trim_galore \
        --length {params.min_length} \
        --output "{params.outdir}" \
        --retain_unpaired \
        --paired \
        -j {threads} \
        -r1 {params.min_length_single} \
        -r2 {params.min_length_single} \
        {input.r1} \
        {input.r2} 2> >(tee {log.errfile} >&2)
        """

