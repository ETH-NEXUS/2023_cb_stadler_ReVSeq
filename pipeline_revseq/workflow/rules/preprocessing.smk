## Helper functions
import os
import pandas as pd

def find_input_fastq_lanes_r1(wildcards):
    #"return a list of all input fastq files available for one sample. Handles multi-lane data"
    mate = "R1"
    sample = str(wildcards.sample)
    matches = pd.read_table(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv", dtype=str)
    snumber = (matches[matches["ethid"] == sample]["Sample number"].to_string(index=False).strip())
    sample_table = sample_map[sample_map["sample"] == sample]
    sample_lanes = sample_table["lane"].tolist()
    sample_lanes.sort()
    fastq_files = []
    total_samples = os.listdir(config["inputOutput"]["input_fastqs"]+"/"+config["plate"])
    for lane in sample_lanes:
        all_samples = [ file for file in total_samples if str(snumber)+"_" in file ]
        lanefile = [ file for file in all_samples if lane in file ]
        raw = [ file for file in lanefile if mate in file ]
        if len(raw) == 0:
            sys.exit("ERROR: can't find raw file for sample " + str(snumber) + " lane " + lane + " mate " + mate)
        fastq_files.append(config["inputOutput"]["input_fastqs"]+"/"+config["plate"]+"/"+raw[0])
    return fastq_files


def find_input_fastq_lanes_r2(wildcards):
    #"return a list of all input fastq files available for one sample. Handles multi-lane data"
    mate = "R2"
    sample = str(wildcards.sample)
    matches = pd.read_table(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv", dtype=str)
    snumber = (matches[matches["ethid"] == sample]["Sample number"].to_string(index=False))
    sample_table = sample_map[sample_map["sample"] == sample]
    sample_lanes = sample_table["lane"].tolist()
    sample_lanes.sort()
    fastq_files = []
    total_samples = os.listdir(config["inputOutput"]["input_fastqs"]+"/"+config["plate"])
    for lane in sample_lanes:
        all_samples = [ file for file in total_samples if str(snumber)+"_" in file ]
        lanefile = [ file for file in all_samples if lane in file ]
        raw = [ file for file in lanefile if mate in file ]
        if len(raw) == 0:
            sys.exit("ERROR: can't find raw file for sample " + str(snumber) + " lane " + lane + " mate " + mate)
        fastq_files.append(config["inputOutput"]["input_fastqs"]+"/"+config["plate"]+"/"+raw[0])
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
    conda:
        "../envs/python.yaml"
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


rule cutadapt:
    input:
        r1 = rules.merge_lanes.output.r1,
        r2 = rules.merge_lanes.output.r2,
    output:
        r1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/cutadapt/{sample}_R1_trimmed.fq.gz"),#temp
        r2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/cutadapt/{sample}_R2_trimmed.fq.gz"),#temp
    params:
        min_length = config["tools"]["cutadapt"]["min_length"],
        adapter_default = config["tools"]["cutadapt"]["adapter_default"],
        samplesheet_dir = config["resources"]["metadata_dir"]+"/"+config["plate"],
        samplesheet_search_string = config["tools"]["cutadapt"]["samplesheet_search_string"],
        anonymization_table = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv",
        sample = "{sample}",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cutadapt/{sample}_cutadapt.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cutadapt/{sample}_cutadapt.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/cutadapt/{sample}.benchmark"
    conda:
        "../envs/trim_galore.yaml"
    threads: config['tools']['cutadapt']['threads']
    shell:
        """
        if (( $(ls -d {params.samplesheet_dir}/{params.samplesheet_search_string}* | wc -l ) == 1 )); then
            samplesheet=$(ls -d {params.samplesheet_dir}/{params.samplesheet_search_string}*)
        else
            echo "ERROR: found 0 or multiple SampleSheet files!"
            exit 1
        fi
        original_sample_name=$(grep {params.sample} {params.anonymization_table} | awk '{{print $1}}')
        adapter1=$(grep ${{original_sample_name}} "${{samplesheet}}" | awk -F ',' '{{print $5}}')
        adapter2=$(grep ${{original_sample_name}} "${{samplesheet}}" | awk -F ',' '{{print $7}}')
        
        # Sometimes we don't receive adapters for a sample
        if [ ${{adapter1}} ] || [[ ${{adapter2}} ]]; then
            adapter_opt=""
            if [ ${{adapter1}} ]; then
                adapter_opt="${{adapter_opt}} -b ${{adapter1}}"
            fi
            if [ ${{adapter2}} ]; then
                adapter_opt="${{adapter_opt}} -b ${{adapter2}}"
            fi
            cutadapt \
            ${{adapter_opt}} -b {params.adapter_default} \
            -O {params.min_length} \
            -o {output.r1} \
            -p {output.r2} \
            -j {threads} \
            -m {params.min_length} \
            {input.r1} \
            {input.r2} 2> >(tee {log.errfile} >&2)
        else
            cutadapt \
            -b {params.adapter_default} \
            -O {params.min_length} \
            -o {output.r1} \
            -p {output.r2} \
            -j {threads} \
            -m {params.min_length} \
            {input.r1} \
            {input.r2} 2> >(tee {log.errfile} >&2)
        fi
        """

#rule trim_primers:
#    input:
#        r1 = rules.merge_lanes.output.r1,
#        r2 = rules.merge_lanes.output.r2,
#    output:
#        r1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_primers/{sample}_merged_R1_val_1.fq.gz"),#temp
#        r2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_primers/{sample}_merged_R2_val_2.fq.gz"),#temp
#    params:
#        base_name = "data/{sample}",
#        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_primers",
#        min_length = config["tools"]["trim_galore"]["min_length"],
#        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
#    log:
#        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/trim_primers/{sample}_trim_galore.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/trim_primers/{sample}_trim_galore.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/trim_primers/{sample}.benchmark"
#    conda:
#        "../envs/trim_galore.yaml"
#    threads: config['tools']['trim_galore']['threads']
#    shell:
#        """
#        trim_galore \
#        --length {params.min_length} \
#        --output "{params.outdir}" \
#        --retain_unpaired \
#        --paired \
#        -j {threads} \
#        -r1 {params.min_length_single} \
#        -r2 {params.min_length_single} \
#        {input.r1} \
#        {input.r2} 2> >(tee {log.errfile} >&2)
#        """
#
#rule trim_adapters:
#    input:
#        r1 = rules.trim_primers.output.r1,
#        r2 = rules.trim_primers.output.r2,
#    output:
#        r1_preliminary = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_adapters/{sample}_merged_R1_val_1_val_1.fq.gz"),#temp
#        r2_preliminary = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_adapters/{sample}_merged_R2_val_2_val_2.fq.gz"),#temp
#        r1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_adapters/{sample}_merged_R1_val_1_val_1_val_1.fq.gz"),#temp
#        r2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_adapters/{sample}_merged_R2_val_2_val_2_val_2.fq.gz"),#temp
#    params:
#        base_name = "data/{sample}",
#        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/trim_adapters",
#        min_length = config["tools"]["trim_galore"]["min_length"],
#        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
#        adapter1 = XXX,
#        adapter2 = XXX,
#    log:
#        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/trim_adapters/{sample}_trim_galore.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/trim_adapters/{sample}_trim_galore.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/trim_adapters/{sample}.benchmark"
#    conda:
#        "../envs/trim_galore.yaml"
#    threads: config['tools']['trim_galore']['threads']
#    shell:
#        """
#        trim_galore \
#        --length {params.min_length} \
#        --output "{params.outdir}" \
#        --retain_unpaired \
#        --paired \
#        -j {threads} \
#        -r1 {params.min_length_single} \
#        -r2 {params.min_length_single} \
#        -a {params.adapter1} \
#        {input.r1} \
#        {input.r2} 2> >(tee {log.errfile} >&2)
#
#        trim_galore \
#        --length {params.min_length} \
#        --output "{params.outdir}" \
#        --retain_unpaired \
#        --paired \
#        -j {threads} \
#        -r1 {params.min_length_single} \
#        -r2 {params.min_length_single} \
#        -a {params.adapter2} \
#        {output.r1_preliminary} \
#        {input.r2_preliminary} 2> >(tee -a {log.errfile} >&2)
#        """