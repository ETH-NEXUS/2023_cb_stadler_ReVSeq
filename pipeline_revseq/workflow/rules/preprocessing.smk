def find_input_fastq_lanes_r1(wildcards):
    "return a list of all input fastq files available for one sample. Handles multi-lane data"
    sample = wildcards.sample
    sample_table = sample_map[sample_map["sample"] == sample]
    sample_lanes = sample_table["lane"].tolist()
    fastq_files_r1 = []
    for lane in sample_lanes:
        fastq_files_r1.append(config["inputOutput"]["input_fastqs"]+ "/" + sample + "/" + sample + "_" + lane + "_R1_001.fastq.gz")
    return fastq_files_r1

def find_input_fastq_lanes_r2(wildcards):
    "return a list of all input fastq files available for one sample. Handles multi-lane data"
    sample = wildcards.sample
    sample_table = sample_map[sample_map["sample"] == sample]
    sample_lanes = sample_table["lane"].tolist()
    fastq_files_r2 = []
    for lane in sample_lanes:
        fastq_files_r2.append(config["inputOutput"]["input_fastqs"]+ "/" + sample + "/" + sample + "_" + lane + "_R2_001.fastq.gz")
    return fastq_files_r2

rule merge_lanes:
    input:
        r1 = find_input_fastq_lanes_r1,
        r2 = find_input_fastq_lanes_r2,
        #r1 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_{lane}_R1_001.fastq.gz", zip, lane=lane_ids, sample = sample_ids),
        #r2 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_{lane}_R2_001.fastq.gz", zip, lane=lane_ids, sample = sample_ids)
    output:
        r1 = ("config["inputOutput"]["output_dir"]/{sample}/merge_lanes/{sample}_merged_R1.fastq.gz"),
        r2 = ("config["inputOutput"]["output_dir"]/{sample}/merge_lanes/{sample}_merged_R2.fastq.gz"),
    log:
        outfile="config["inputOutput"]["output_dir"]/logs/{sample}/merge_lanes/merge_lanes.out.log",
        errfile="config["inputOutput"]["output_dir"]/logs/{sample}/merge_lanes/merge_lanes.err.log",
    benchmark:
        "config["inputOutput"]["output_dir"]/logs/benchmark/{sample}/merge_lanes/merge_lanes.benchmark"
    shell:
        "gzcat {input.r1} | gzip > {output.r1} && gzcat {input.r2} | gzip > {output.r2}"


#rule merge_refs:
#    input:
#        refs = config["resources"]["reference"],
#        host_ref = config["resources"]["host_ref"],
#    output:
#        outref = "config["inputOutput"]["output_dir"]/merge_refs/ref.fa",
#    log:
#        outfile="config["inputOutput"]["output_dir"]/logs/merge_refs/merge_refs.out.log",
#        errfile="config["inputOutput"]["output_dir"]/logs/merge_refs/merge_refs.err.log",
#    benchmark:
#        "config["inputOutput"]["output_dir"]/logs/benchmark/merge_refs/merge_refs.benchmark"
#    conda:
#        "../envs/consensus.yaml"
#    shell:
#        "python workflow/scripts/merge_refs.py --refs {input.refs} --host_ref {input.host_ref} --output {output.outref}"


rule bwa_index:
    input:
        reference = config["resources"]["reference"]
        #reference = rules.merge_refs.output.outref
    output:
        referenceout = config["resources"]["reference"] + '.bwt'
        #referenceout = rules.merge_refs.output.outref + '.bwt'
    log:
        outfile="config["inputOutput"]["output_dir"]/logs/bwa_index/bwa_index.out.log",
        errfile="config["inputOutput"]["output_dir"]/logs/bwa_index/bwa_index.err.log",
    benchmark:
        "config["inputOutput"]["output_dir"]/logs/benchmark/bwa_index/bwa_index.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.reference}"


rule trim_galore:
    input:
        r1 = rules.merge_lanes.output.r1,
        r2 = rules.merge_lanes.output.r2,
    output:
        r1 = "config["inputOutput"]["output_dir"]/{sample}/trim_galore/{sample}_merged_R1_val_1.fq.gz",
        r2 = "config["inputOutput"]["output_dir"]/{sample}/trim_galore/{sample}_merged_R2_val_2.fq.gz",
    params:
        base_name = "data/{sample}",
        outdir = 'config["inputOutput"]["output_dir"]/{sample}/trim_galore',
        min_length = config["tools"]["trim_galore"]["min_length"],
        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
    log:
        outfile="config["inputOutput"]["output_dir"]/logs/{sample}/trim_galore/{sample}_trim_galore.out.log",
        errfile="config["inputOutput"]["output_dir"]/logs/{sample}/trim_galore/{sample}_trim_galore.err.log",
    benchmark:
        "config["inputOutput"]["output_dir"]/logs/benchmark/trim_galore/{sample}.benchmark"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        'trim_galore --length {params.min_length} --output "{params.outdir}" --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {input.r1} {input.r2}'

