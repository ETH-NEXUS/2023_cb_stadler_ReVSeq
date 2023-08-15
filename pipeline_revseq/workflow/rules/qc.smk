rule fastqc_raw:
    input:
        inputdir = expand(config["inputOutput"]["input_fastqs"]+"/"+config["plate"]+"/{{sample}}/")
    output:
        zip1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_raw/{sample}_L001_R1_001_fastqc.zip",
        zip2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_raw/{sample}_L001_R2_001_fastqc.zip"
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_raw"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/fastqc_raw/fastqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/fastqc_raw/fastqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/fastqc_raw/fastqc.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        "fastqc {input.inputdir}/{wildcards.sample}_L00*_R*_001.fastq.gz -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)"


rule fastqc_merged:
    input:
        r1 = rules.merge_lanes.output.r1,
        r2 = rules.merge_lanes.output.r2,
    output:
        zip1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_merged/{sample}_merged_R1_fastqc.zip",
        zip2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_merged/{sample}_merged_R2_fastqc.zip",
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_merged"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/fastqc/fastqc_merged.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/fastqc/fastqc_merged.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/fastqc/fastqc_merged.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        """
        fastqc {input.r1} -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)
        fastqc {input.r2} -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)
        """


rule samtoolsstats:
    input:
        bam = rules.remove_duplicates.output.bam,
    output:
        stats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/samtoolsstats/{sample}_samtools_stats.txt",
        flagstats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/samtoolsstats/{sample}_samtools_flagstats.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bamtools/bamtools.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bamtools/bamtools.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/bamtools/bamtools.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.stats} && samtools flagstats {input.bam} > {output.flagstats} 2> >(tee {log.errfile} >&2)"


rule rseqc:
    input:
        bam = rules.remove_duplicates.output.bam,
        index = rules.samtools_index.output.index
    output:
        stats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/rseqc/{sample}_bam_stat.txt",
    params:
        clipprefix = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/rseqc/{sample}"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/rseqc/rseqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/rseqc/rseqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/rseqc/rseqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output.stats} 2> >(tee {log.errfile} >&2)"


rule multiqc:
    input:
        fastqcresult1 = expand(rules.fastqc_merged.output.zip1, sample=sample_ids),
        fastqcresult2 = expand(rules.fastqc_merged.output.zip2, sample=sample_ids),
        trimresult1 =  expand(rules.trim_galore.output.r1, sample=sample_ids),
        trimresult2 =  expand(rules.trim_galore.output.r2, sample=sample_ids),
        stats = expand(rules.samtoolsstats.output.stats, sample=sample_ids),
        flagstats = expand(rules.samtoolsstats.output.flagstats, sample=sample_ids),
        rseqcstats = expand(rules.rseqc.output.stats, sample=sample_ids),
        qualimap = expand(rules.qualimap.output.report, sample=sample_ids),
    output:
        outfile = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc/multiqc_report.html"
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc",
        inputdir = config["inputOutput"]["output_dir"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc/multiqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc/multiqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/multiqc/multiqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.inputdir} --ignore "{params.inputdir}/*/fastqc_raw" --ignore "{params.inputdir}/*/dh_fastqc" -o {params.outdir}  2> >(tee {log.errfile} >&2)
        """
