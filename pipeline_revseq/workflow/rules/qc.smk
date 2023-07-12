rule fastqc:
    input:
        inputdir = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/")
    #    r1 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_L00*_R1_001.fastq.gz", lane=lane_ids, sample = sample_ids),
    #    r2 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_L00*_R2_001.fastq.gz", lane=lane_ids, sample = sample_ids)
    output:
        zip1 = "results/{sample}/fastqc/{sample}_L001_R1_001_fastqc.zip",
        zip2 = "results/{sample}/fastqc/{sample}_L001_R2_001_fastqc.zip"
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
        #"fastqc {input.r1} {input.r2} -t {threads} -o {params.outdir}"
        "fastqc {input.inputdir}/{wildcards.sample}_L00*_R*_001.fastq.gz -t {threads} -o {params.outdir}"


rule samtoolsstats:
    input:
        bam = rules.filter_host_reads.output.bam
    output:
        stats = "results/{sample}/samtoolsstats/samtools_stats.txt",
        flagstats = "results/{sample}/samtoolsstats/samtools_flagstats.txt",
    log:
        outfile="logs/{sample}/bamtools/bamtools.out.log",
        errfile="logs/{sample}/bamtools/bamtools.err.log",
    benchmark:
        "logs/benchmark/{sample}/bamtools/bamtools.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.stats} && samtools flagstats {input.bam} > {output.flagstats}"


rule rseqc:
    input:
        bam = "results/{sample}/bwa/mapped_reads.bam"
    output:
        stats = "results/{sample}/rseqc/bam_stat.txt",
        clipr1 = "results/{sample}/rseqc/{sample}.clip_profile.R1.pdf",
        clipr2 = "results/{sample}/rseqc/{sample}.clip_profile.R2.pdf",
    params:
        clipprefix = "results/{sample}/rseqc/{sample}"
    log:
        outfile="logs/{sample}/rseqc/rseqc.out.log",
        errfile="logs/{sample}/rseqc/rseqc.err.log",
    benchmark:
        "logs/benchmark/{sample}/rseqc/rseqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        bam_stat.py -i {input.bam} > {output.stats}
        clipping_profile.py -i {input.bam} -s "PE" -o {params.clipprefix}
        """


rule multiqc:
    input:
        fastqcresult1 = rules.fastqc.output.zip1,
        fastqcresult2 = rules.fastqc.output.zip2,
        trimresult1 =  rules.trim_galore.output.r1,
        trimresult2 =  rules.trim_galore.output.r2,
        stats = rules.samtoolsstats.output.stats,
        flagstats = rules.samtoolsstats.output.flagstats,
        rseqcstats = rules.rseqc.output.stats,
        rseqcclip1 = rules.rseqc.output.clipr1,
        rseqcclip2 = rules.rseqc.output.clipr2,
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
