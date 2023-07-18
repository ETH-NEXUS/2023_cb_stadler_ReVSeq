rule fastqc:
    input:
        inputdir = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/")
    #    r1 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_L00*_R1_001.fastq.gz", lane=lane_ids, sample = sample_ids),
    #    r2 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_L00*_R2_001.fastq.gz", lane=lane_ids, sample = sample_ids)
    output:
        zip1 = config["inputOutput"]["output_dir"]+"/{sample}/fastqc/{sample}_L001_R1_001_fastqc.zip",
        zip2 = config["inputOutput"]["output_dir"]+"/{sample}/fastqc/{sample}_L001_R2_001_fastqc.zip"
    params:
        outdir = config["inputOutput"]["output_dir"]+"/{sample}/fastqc"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/fastqc/fastqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/fastqc/fastqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/{sample}/fastqc/fastqc.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        #"fastqc {input.r1} {input.r2} -t {threads} -o {params.outdir}"
        "fastqc {input.inputdir}/{wildcards.sample}_L00*_R*_001.fastq.gz -t {threads} -o {params.outdir}"


rule samtoolsstats:
    input:
        #bam = rules.filter_host_reads.output.bam
        bam = rules.sort.output.outfile
    output:
        stats = config["inputOutput"]["output_dir"]+"/{sample}/samtoolsstats/{sample}_samtools_stats.txt",
        flagstats = config["inputOutput"]["output_dir"]+"/{sample}/samtoolsstats/{sample}_samtools_flagstats.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/bamtools/bamtools.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/bamtools/bamtools.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/{sample}/bamtools/bamtools.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.stats} && samtools flagstats {input.bam} > {output.flagstats}"


rule rseqc:
    input:
        bam = rules.sort.output.outfile,
        index = rules.samtools_index.output.index
    output:
        stats = config["inputOutput"]["output_dir"]+"/{sample}/rseqc/{sample}_bam_stat.txt",
    params:
        clipprefix = config["inputOutput"]["output_dir"]+"/{sample}/rseqc/{sample}"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/rseqc/rseqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/rseqc/rseqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/{sample}/rseqc/rseqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output.stats}"


rule multiqc:
    input:
        fastqcresult1 = expand(rules.fastqc.output.zip1, sample=sample_ids),
        fastqcresult2 = expand(rules.fastqc.output.zip2, sample=sample_ids),
        trimresult1 =  expand(rules.trim_galore.output.r1, sample=sample_ids),
        trimresult2 =  expand(rules.trim_galore.output.r2, sample=sample_ids),
        stats = expand(rules.samtoolsstats.output.stats, sample=sample_ids),
        flagstats = expand(rules.samtoolsstats.output.flagstats, sample=sample_ids),
        rseqcstats = expand(rules.rseqc.output.stats, sample=sample_ids),
    output:
        outfile = config["inputOutput"]["output_dir"]+"/multiqc/multiqc_report.html"
    params:
        outdir = config["inputOutput"]["output_dir"]+"/multiqc",
        inputdir = config["inputOutput"]["output_dir"]
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/multiqc/multiqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/multiqc/multiqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/multiqc/multiqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc {params.inputdir} -o {params.outdir}"
