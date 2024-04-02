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
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/fastqc/{sample}_fastqc_merged.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        """
        fastqc {input.r1} -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)
        fastqc {input.r2} -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)
        """


rule fastqc_trimmed:
    input:
        r1 = rules.cutadapt.output.r1,
        r2 = rules.cutadapt.output.r2,
    output:
        zip1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_trimmed/{sample}_R1_trimmed_fastqc.zip",
        zip2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_trimmed/{sample}_R2_trimmed_fastqc.zip",
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/fastqc_trimmed"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/fastqc/fastqc_trimmed.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/fastqc/fastqc_trimmed.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/fastqc/{sample}fastqc_trimmed.benchmark"
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
        #bam = rules.bwa.output.bam,
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam"
    output:
        stats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/samtoolsstats/{sample}_samtools_stats.txt",
        flagstats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/samtoolsstats/{sample}_samtools_flagstats.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/samtoolsstats/samtoolsstats.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/samtoolsstats/samtoolsstats.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/samtoolsstats/{sample}_samtoolsstats.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.stats} && samtools flagstats {input.bam} > {output.flagstats} 2> >(tee {log.errfile} >&2)"


rule rseqc:
    input:
        #bam = rules.bwa.output.bam,
        #index = rules.bwa.output.bai
        bam =config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam",
        index = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam.bai"
    output:
        stats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/rseqc/{sample}_bam_stat.txt",
    params:
        clipprefix = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/rseqc/{sample}"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/rseqc/rseqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/rseqc/rseqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/rseqc/{sample}_rseqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output.stats} 2> >(tee {log.errfile} >&2)"


rule qualimap:
    input:
        #bam = rules.bwa.output.bam,
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam"
    output:
        report = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap/qualimapReport.html",
        genome_res = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap/genome_results.txt",
    params:
        regions = config["resources"]["reference_table"],
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap/",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/qualimap/qualimap.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/qualimap/qualimap.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/qualimap/{sample}_qualimap.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        unset DISPLAY
        qualimap bamqc -outdir {params.outdir} -bam {input.bam} --feature-file {params.regions} -c 2> >(tee {log.errfile} >&2)
        """


rule multiqc:
    input:
        fastqcresult1 = expand(rules.fastqc_merged.output.zip1, sample=sample_ids),
        fastqcresult2 = expand(rules.fastqc_merged.output.zip2, sample=sample_ids),
        stats = expand(rules.samtoolsstats.output.stats, sample=sample_ids),
        flagstats = expand(rules.samtoolsstats.output.flagstats, sample=sample_ids),
        rseqcstats = expand(rules.rseqc.output.stats, sample=sample_ids),
    output:
        outfile = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc/multiqc_report.html",
        outdir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc")
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc",
        inputdir = config["inputOutput"]["output_dir"]+"/"+config["plate"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc/multiqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc/multiqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/multiqc/multiqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.inputdir} \
        --ignore "{params.inputdir}/*/fastqc_raw" \
        --ignore "{params.inputdir}/*/fastqc_trimmed" \
        --ignore "{params.inputdir}/*/dh_fastqc" \
        --ignore "{params.inputdir}/*/qualimap_filtered" \
        --ignore "{params.inputdir}/*/rseqc_filtered" \
        --ignore "{params.inputdir}/*/samtoolsstats_filtered" \
        --ignore "{params.inputdir}/*/qualimap" \
        --interactive \
        -o {params.outdir}  2> >(tee {log.errfile} >&2)
        """


rule multiqc_trimmed:
    input:
        trimresult1 =  expand(rules.fastqc_trimmed.output.zip1, sample=sample_ids),
        trimresult2 =  expand(rules.fastqc_trimmed.output.zip2, sample=sample_ids),
        stats = expand(rules.samtoolsstats.output.stats, sample=sample_ids),
        flagstats = expand(rules.samtoolsstats.output.flagstats, sample=sample_ids),
        rseqcstats = expand(rules.rseqc.output.stats, sample=sample_ids),
    output:
        outfile = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc_trimmed/multiqc_report.html",
        outdir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc_trimmed")
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc_trimmed",
        inputdir = config["inputOutput"]["output_dir"]+"/"+config["plate"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc_trimmed/multiqc_trimmed.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc_trimmed/multiqc_trimmed.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/multiqc_trimmed/multiqc_trimmed.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.inputdir} \
        --ignore "{params.inputdir}/*/fastqc_raw" \
        --ignore "{params.inputdir}/*/fastqc_merged" \
        --ignore "{params.inputdir}/*/dh_fastqc" \
        --ignore "{params.inputdir}/*/qualimap_filtered" \
        --ignore "{params.inputdir}/*/rseqc_filtered" \
        --ignore "{params.inputdir}/*/samtoolsstats_filtered" \
        --ignore "{params.inputdir}/*/qualimap" \
        --interactive \
        -o {params.outdir}  2> >(tee {log.errfile} >&2)
        """


rule samtoolsstats_filtered:
    input:
        bam = rules.remove_duplicates.output.bam,
    output:
        stats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/samtoolsstats_filtered/{sample}_samtools_stats.txt",
        flagstats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/samtoolsstats_filtered/{sample}_samtools_flagstats.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/samtoolsstats_filtered/samtoolsstats.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/samtoolsstats_filtered/samtoolsstats.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/samtoolsstats_filtered/{sample}_samtoolsstats.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.stats} && samtools flagstats {input.bam} > {output.flagstats} 2> >(tee {log.errfile} >&2)"


rule rseqc_filtered:
    input:
        bam = rules.remove_duplicates.output.bam,
        index = rules.samtools_index.output.index
    output:
        stats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/rseqc_filtered/{sample}_bam_stat.txt",
    params:
        clipprefix = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/rseqc_filtered/{sample}"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/rseqc_filtered/rseqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/rseqc_filtered/rseqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/rseqc_filtered/{sample}_rseqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        "bam_stat.py -i {input.bam} > {output.stats} 2> >(tee {log.errfile} >&2)"


rule multiqc_filtered:
    input:
        stats = expand(rules.samtoolsstats_filtered.output.stats, sample=sample_ids),
        flagstats = expand(rules.samtoolsstats_filtered.output.flagstats, sample=sample_ids),
        rseqcstats = expand(rules.rseqc_filtered.output.stats, sample=sample_ids),
        #qualimap = expand(rules.qualimap_filtered.output.report, sample=sample_ids),
    output:
        outfile = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc_filtered/multiqc_report.html",
        outdir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc_filtered")
    params:
        inputdir = config["inputOutput"]["output_dir"]+"/"+config["plate"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc_filtered/multiqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/multiqc_filtered/multiqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/multiqc_filtered/multiqc.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc {params.inputdir} \
        --ignore "{params.inputdir}/*/fastqc_raw" \
        --ignore "{params.inputdir}/*/dh_fastqc" \
        --ignore "{params.inputdir}/*/fastqc_merged" \
        --ignore "{params.inputdir}/*/cutadapt" \
        --ignore "{params.inputdir}/*/qualimap" \
        --ignore "{params.inputdir}/*/qualimap_filtered" \
        --ignore "{params.inputdir}/*/rseqc" \
        --ignore "{params.inputdir}/*/samtoolsstats" \
        --ignore "{params.inputdir}/*/remove_duplicates" \
        --interactive \
        -o {output.outdir}  2> >(tee {log.errfile} >&2)
        """