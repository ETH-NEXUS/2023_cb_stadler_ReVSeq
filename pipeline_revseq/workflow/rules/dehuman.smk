rule bwa:
    input:
        ref = rules.merge_refs.output.referenceout,
        ref_index = rules.bwa_index.output.ref_index,
        r1 = rules.cutadapt.output.r1,
        r2 = rules.cutadapt.output.r2,
    output:
        bwa = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_mapped_reads.bam"),
        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam"),
        bai = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam.bai"),
        readcount_all = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_readcount_all.txt",
        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_readcount.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bwa/bwa.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bwa/bwa.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/bwa/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    threads: config["tools"]["bwa"]["threads"]
    resources:
        mem_mb=6500
    shell:
        """
        bwa mem -t {threads} -M {input.ref} {input.r1} {input.r2} > {output.bwa} 2> >(tee {log.errfile} >&2)
        samtools sort -@ {threads} --output-fmt=SAM {output.bwa} | samtools view -Sb - > {output.bam} 2> >(tee -a {log.errfile} >&2)
        samtools index {output.bam} 2> >(tee -a {log.errfile} >&2)
        samtools view -c {output.bam} > {output.readcount_all} 2> >(tee -a {log.errfile} >&2)
        samtools view -c -F 4 {output.bam} > {output.readcount} 2> >(tee -a {log.errfile} >&2)
        """


rule dehuman:
    input:
        bam = rules.bwa.output.bam,
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam"),
        bai = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam.bai"),
        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_readcount.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dehuman/dehuman.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dehuman/dehuman.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/dehuman/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -h {input.bam} | grep -v -e "chr" | samtools view -Sbh > {output.bam} 
        samtools index {output.bam}
        samtools view -c -F 4 {output.bam} > {output.readcount}
        """


rule cram:
    input:
        bam = rules.dehuman.output.bam,
        ref = rules.merge_refs.output.referenceout,
    output:
        cram = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/cram/{sample}.cram",
    params:
	    opt = config["tools"]["cram"]["tool_options"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cram/fastq_to_cram.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cram/fastq_to_cram.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/cram/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools sort -T {input.bam} --reference {input.ref} --output-fmt {params.opt} -o {output.cram} {input.bam}"


rule flat_file:
    input:
        cram = rules.cram.output.cram,
    output:
        fasta = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/flat_file/{sample}.fasta",
        flat_temp = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/flat_file/{sample}.embl.tmp"),
        flat = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/flat_file/{sample}.embl",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/flat_file/flat_file.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/flat_file/flat_file.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/flat_file/{sample}.benchmark"
    conda:
        "../envs/seqret.yaml"
    shell:
    """
        samtools fasta {input.cram} > {output.fasta}
        seqret -sequence {output.fasta} -outseq {output.flat_temp} -osformat embl
        #The flat file does not fully comply with the upload requirements. The header needs to be manually modified
        cp {output.flat_temp} {output.flat}
    """

rule bam_to_fastq:
    input:
        bam = rules.dehuman.output.bam,
    output:
        fq1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/{sample}_bam_to_fastq_r1.fastq",#temp
        fq2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/{sample}_bam_to_fastq_r2.fastq",#temp
        outdir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/"),
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bam_to_fastq/bam_to_fastq.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bam_to_fastq/bam_to_fastq.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/bam_to_fastq/{sample}_bam_to_fastq.benchmark"
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
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/dh_fastqc/{sample}_dh_fastqc.benchmark"
    conda:
        "../envs/qc.yaml"
    threads: config["threads"]
    shell:
        """
        fastqc {input.inputdir}/{wildcards.sample}_bam_to_fastq_r1.fastq -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)
        echo test
        fastqc {input.inputdir}/{wildcards.sample}_bam_to_fastq_r2.fastq -t {threads} -o {params.outdir}  2> >(tee {log.errfile} >&2)
        """

