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
        r1 = "results/{sample}/merge_lanes/merged_R1.fastq.gz",
        r2 = "results/{sample}/merge_lanes/merged_R2.fastq.gz",
    log:
        outfile="logs/{sample}/merge_lanes/merge_lanes.out.log",
        errfile="logs/{sample}/merge_lanes/merge_lanes.err.log",
    benchmark:
        "logs/benchmark/{sample}/merge_lanes/merge_lanes.benchmark"
    shell:
        "gzcat {input.r1} | gzip > {output.r1} && gzcat {input.r2} | gzip > {output.r2}"


#rule merge_refs:
#    input:
#        refs = config["resources"]["reference"],
#        host_ref = config["resources"]["host_ref"],
#    output:
#        outref = "results/merge_refs/ref.fa",
#    log:
#        outfile="logs/merge_refs/merge_refs.out.log",
#        errfile="logs/merge_refs/merge_refs.err.log",
#    benchmark:
#        "logs/benchmark/merge_refs/merge_refs.benchmark"
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
        outfile="logs/bwa_index/bwa_index.out.log",
        errfile="logs/bwa_index/bwa_index.err.log",
    benchmark:
        "logs/benchmark/bwa_index/bwa_index.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.reference}"


rule trim_galore:
    input:
        r1 = rules.merge_lanes.output.r1,
        r2 = rules.merge_lanes.output.r2,
    output:
        r1 = "results/{sample}/trim_galore/merged_R1_val_1.fq.gz",
        r2 = "results/{sample}/trim_galore/merged_R2_val_2.fq.gz",
    params:
        base_name = "data/{sample}",
        outdir = 'results/{sample}/trim_galore',
        min_length = config["tools"]["trim_galore"]["min_length"],
        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
    log:
        outfile="logs/{sample}/trim_galore/trim_galore.out.log",
        errfile="logs/{sample}/trim_galore/trim_galore.err.log",
    benchmark:
        "logs/benchmark/trim_galore/{sample}.benchmark"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        'trim_galore --length {params.min_length} --output "{params.outdir}" --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {input.r1} {input.r2}'

