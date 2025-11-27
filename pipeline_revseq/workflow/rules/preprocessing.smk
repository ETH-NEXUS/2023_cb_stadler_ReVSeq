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
        r1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_lanes/{sample}_merged_R1.fastq.gz",
        r2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_lanes/{sample}_merged_R2.fastq.gz",
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

rule kraken2_cp:
    input:
        kraken2_report = f"/data/raw_data/results/{config['plate']}/{{sample}}/{{sample}}.report"
    output:
        kraken2_report_cp = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/kraken2/{sample}.report",
    shell:
        "cp {input.kraken2_report} {output.kraken2_report_cp}"


rule kraken2:
    input:
        r1 =rules.merge_lanes.output.r1,
        r2 =rules.merge_lanes.output.r2,
    output:
        kraken2_report = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/kraken2/{sample}.report",
    params:
        db = config["resources"]["kraken2_db"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/kraken2/kraken2.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/kraken2/kraken2.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/kraken2/{sample}_kraken2.benchmark"
    conda:
        "../envs/kraken2.yaml"
    threads: config['tools']['kraken2']['threads']
    shell:
        """
        kraken2 --db {params.db} \
            --paired \
            --report {output.kraken2_report} \
            --use-names \
            --threads {threads} \
            {input.r1} {input.r2} 2> >(tee {log.errfile} >&2)
        """

rule kraken2_report:
    input:       
        kraken2_report = rules.kraken2.output.kraken2_report,
    output:
        detected_substrains = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/kraken2/{sample}.detected_substrains.tsv",
    params:
        k2 = "/data/k2/"+config["plate"]+"/{sample}/{sample}.report",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/kraken2/kraken2.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/kraken2/kraken2.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/kraken2/{sample}_kraken2.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        
        python workflow/scripts/kraken_report_parser.py \
            --report {input.kraken2_report} \
            --output {output.detected_substrains}
        """


rule gather_references:
    input:
        detected_substrains = rules.kraken2.output.detected_substrains,
        missing=config["resources"]["reference_dir"]+"/missing.txt ",
    output:
        bed = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/gather_references/{sample}.bed",
        fasta_links_dir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/gather_references"),
    params:
        fastadir = config["resources"]["reference_dir"]+"/fasta",
        download_list = config["resources"]["reference_dir"]+"/library_report.tsv",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/gather_references/gather_references.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/gather_references/gather_references.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/gather_references/{sample}_gather_references.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/gather_references.py \
            --substrains {input.detected_substrains} \
            --downloadlist {params.download_list} \
            --fastadir {params.fastadir} \
            --bed {output.bed} \
            --missing {input.missing} \
            --fastalinksdir {output.fasta_links_dir}
        """


rule merge_refs:
    input:
        bed = rules.gather_references.output.bed,
    output:
        referenceout = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_refs/{sample}_merged_virus_host_ref.fa",
        referenceout_virus_only = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/merge_refs/{sample}_merged_virus_only_ref.fa",
    params:
        fasta_links_dir = rules.gather_references.output.fasta_links_dir,
        host_virus_reference = config["resources"]["host_ref"],
        virus_reference = config["resources"]["reference"]
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/merge_refs/merge_refs.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/merge_refs/merge_refs.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/merge_refs/{sample}_merge_refs.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        if [[ -f {output.referenceout} ]]; then
            rm {output.referenceout}
        fi
        if [[ -f {output.referenceout_virus_only} ]];then
            rm {output.referenceout_virus_only}
        fi
        if find {params.fasta_links_dir} -maxdepth 1 -type f -name "*.gz" | grep -q .; then
            zcat {params.fasta_links_dir}/*gz >> {output.referenceout} 2> >(tee {log.errfile} >&2)
            cp {output.referenceout} {output.referenceout_virus_only} 2> >(tee {log.errfile} >&2)
        fi        
        cat {params.host_virus_reference} >> {output.referenceout} 2> >(tee {log.errfile} >&2)
        cat {params.virus_reference} >> {output.referenceout_virus_only} 2> >(tee {log.errfile} >&2)
        """


rule bwa_index:
    input:
        reference = rules.merge_refs.output.referenceout
    output:
        ref_index = rules.merge_refs.output.referenceout + '.bwt'
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/bwa_index/bwa_index.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/bwa_index/bwa_index.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/bwa_index/{sample}_bwa_index.benchmark"
    resources:
        mem_mb = 5000,
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
        "../envs/cutadapt.yaml"
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
