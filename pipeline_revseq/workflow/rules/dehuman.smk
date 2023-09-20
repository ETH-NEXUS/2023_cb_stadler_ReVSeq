#rule bwa:
#    input:
#        ref = config["resources"]["host_ref"],
#        ref_index = rules.bwa_index.output.ref_index,
#        reads = rules.trim_galore.output
#    output:
#        bwa = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_mapped_reads.bam"),
#        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam"),
#        bai = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam.bai"),
#        readcount_all = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_readcount_all.txt",
#        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_readcount.txt",
#    log:
#        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bwa/bwa.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/bwa/bwa.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/bwa/{sample}.benchmark"
#    conda:
#        "../envs/bwa.yaml"
#    threads: config["threads"]
#    shell:
#        """
#        bwa mem -t {threads} -M {input.ref} {input.reads} > {output.bwa}
#        samtools sort -@ {threads} --output-fmt=SAM {output.bwa} | samtools view -Sb - > {output.bam} 2> >(tee {log.errfile} >&2)
#        samtools index {output.bam}
#        samtools view -c {output.bam} > {output.readcount_all}
#        samtools view -c -F 4 {output.bam} > {output.readcount}
#        """
#
#
#rule dehuman:
#    input:
#        #bam = rules.bwa.output.bam,
#        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}.bam"
#    output:
#        bam = (config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam"),
#        bai = (config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam.bai"),
#        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_readcount.txt",
#    log:
#        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dehuman/dehuman.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/dehuman/dehuman.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/dehuman/{sample}.benchmark"
#    conda:
#        "../envs/samtools.yaml"
#    shell:
#        """
#        samtools view -h {input.bam} | grep -v -e "chr" | samtools view -Sbh > {output.bam} 
#        samtools index {output.bam}
#        samtools view -c -F 4 {output.bam} > {output.readcount}
#        """


rule cram:
    input:
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam",
        #bam = rules.dehuman.output.bam,
    output:
        cram = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/cram/{sample}.cram",
    params:
        ref = config["resources"]["host_ref"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cram/fastq_to_cram.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/cram/fastq_to_cram.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/cram/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -C -T {params.ref} -o {output.cram} {input.bam} 2> >(tee {log.errfile} >&2)"


rule bam_to_fastq:
    input:
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam",
        #bam = rules.dehuman.output.bam,
    output:
        fq1 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/{sample}_bam_to_fastq_r1.fastq"),#temp
        fq2 = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/{sample}_bam_to_fastq_r2.fastq"),#temp
        outdir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bam_to_fastq/"),
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

