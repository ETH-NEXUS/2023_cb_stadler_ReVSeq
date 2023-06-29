rule fastqc:
    input:
        r1 = config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_R1.fastq.gz",
        r2 = config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_R2.fastq.gz",
    output:
        zip1 = "results/{sample}/fastqc/{sample}_R1_fastqc.zip",
        zip2 = "results/{sample}/fastqc/{sample}_R2_fastqc.zip"
    params:
        outdir = "results/{sample}/fastqc"
    log:
        outfile="logs/{sample}/fastqc/fastqc.out.log",
        errfile="logs/{sample}/fastqc/fastqc.err.log",
    benchmark:
        "logs/benchmark/{sample}/fastqc/fastqc.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        "fastqc {input.r1} {input.r2} -t {threads} -o {params.outdir}"


rule samtoolsstats:
    input:
        inputfile = "results/{sample}/bwa/mapped_reads.bam"
    output:
        stats = "results/{sample}/bamqc/samtools_stats.txt",
        flagstats = "results/{sample}/bamqc/samtools_flagstats.txt",
    log:
        outfile="logs/{sample}/bamtools/bamtools.out.log",
        errfile="logs/{sample}/bamtools/bamtools.err.log",
    benchmark:
        "logs/benchmark/{sample}/bamtools/bamtools.benchmark"
    conda:
        "../envs/bwa.yaml"
    threads: config["threads"]
    shell:
        "samtools stats -@ {threads} {input.inputfile} > {output.stats} && samtools flagstats {input.inputfile} > {output.flagstats}"


rule multiqc:
    input:
        fastqcresult1 = expand("results/{sample}/fastqc/{sample}_R1_fastqc.zip", sample=sample_ids),
        fastqcresult2 = expand("results/{sample}/fastqc/{sample}_R2_fastqc.zip", sample=sample_ids),
        trimresult1 =  expand("results/{sample}/trim_galore/{sample}_R1.fastq.gz_trimming_report.txt", sample=sample_ids),
        trimresult2 =  expand("results/{sample}/trim_galore/{sample}_R2.fastq.gz_trimming_report.txt", sample=sample_ids),
        stats = expand("results/{sample}/bamqc/samtools_stats.txt", sample=sample_ids),
        flagstats = expand("results/{sample}/bamqc/samtools_flagstats.txt", sample=sample_ids),
    output:
        outfile = "results/multiqc/multiqc_report.html"
    params:
        outdir = "results/multiqc"
    log:
        outfile="logs/multiqc/multiqc.out.log",
        errfile="logs/multiqc/multiqc.err.log",
    benchmark:
        "logs/benchmark/multiqc/multiqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        "multiqc results -o {params.outdir}"
