rule bwa:
    input:
        ref = rules.merge_refs.output.referenceout,
        index = rules.bwa_index.output.index,
        reads = rules.trim_galore.output
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_mapped_reads.bam")
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bwa/bwa.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bwa/bwa.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/bwa/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem -M {input.ref} {input.reads} |samtools view -Sb - > {output.bam} 2> >(tee {log.errfile} >&2)"


rule dehuman:
    input:
        bam = rules.bwa.output.bam,
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam"),
        only_human = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_only_human.bam"),
    params:
        human_regions = config["tools"]["dehuman"]["human_regions"]
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dehuman/dehuman.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dehuman/dehuman.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/dehuman/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -hb -U {output.bam} -o {output.only_human} {input.bam} {params.human_regions} 2> >(tee {log.errfile} >&2)"


rule cram:
    input:
        bam = rules.dehuman.output.bam,
    output:
        cram = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/cram/{sample}.cram",
    params:
        ref = rules.merge_refs.output.referenceout,
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cram/fastq_to_cram.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cram/fastq_to_cram.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/cram/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -C -T {params.ref} -o {output.cram} {input.bam} 2> >(tee {log.errfile} >&2)"


rule bam_to_fastq:
    input:
        bam = rules.dehuman.output.bam,
    output:
        fq1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/{sample}_bam_to_fastq_r1.fastq"),
        fq2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/{sample}_bam_to_fastq_r2.fastq"),
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bam_to_fastq/bam_to_fastq.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bam_to_fastq/bam_to_fastq.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/bam_to_fastq/bam_to_fastq.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools collate -u -O {input.bam} | samtools fastq -1 {output.fq1} -2 {output.fq2} -0 /dev/null -s /dev/null -n 2> >(tee {log.errfile} >&2)"


rule dh_fastqc:
    input:
        inputdir = rules.bam_to_fastq.output.outdir,
    output:
        zip1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dh_fastqc/{sample}_bam_to_fastq_r1_fastqc.zip",
        zip2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dh_fastqc/{sample}_bam_to_fastq_r2_fastqc.zip"
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dh_fastqc"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dh_fastqc/dh_fastqc.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dh_fastqc/dh_fastqc.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/dh_fastqc/dh_fastqc.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        "fastqc {input.inputdir}/{wildcards.sample}_bam_to_fastq_r*.fastq -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)"

