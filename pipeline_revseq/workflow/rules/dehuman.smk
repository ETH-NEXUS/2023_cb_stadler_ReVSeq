rule fasta_index:
    input:
        reference = config["resources"]["reference"]
    output:
        outindex = config["resources"]["reference"] + ".faidx"
    log:
        outfile="logs/fasta_index/fasta_index.out.log",
        errfile="logs/fasta_index/fasta_index.err.log",
    benchmark:
        "logs/benchmark/fasta_index/fasta_index.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.reference} --fai-idx {output.outindex}"


rule bwa:
    input:
        ref = config["resources"]["reference"],
        #ref = rules.merge_refs.output.outref,
        index = rules.bwa_index.output.referenceout,
        reads = rules.trim_galore.output
    output:
        bam = "results/{sample}/bwa/mapped_reads.bam"
    log:
        outfile="logs/{sample}/bwa/bwa.out.log",
        errfile="logs/{sample}/bwa/bwa.err.log",
    benchmark:
        "logs/benchmark/bwa/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} |samtools view -Sb - > {output.bam}"


rule samtools_index:
    input:
        bam = rules.bwa.output.bam
    output:
        index = rules.bwa.output.bam + '.bai'
    log:
        outfile="logs/{sample}/bwa/samtools_index.out.log",
        errfile="logs/{sample}/bwa/samtools_index.err.log",
    benchmark:
        "logs/benchmark/bwa/{sample}_samtools_index.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bam}"


rule filter_host_reads:
    input:
        bam = rules.bwa.output.bam,
        refindex = rules.fasta_index.output.outindex,
    output:
        bam = temp("results/{sample}/filter_host_reads/{sample}_dehuman.bam"),
        cram = "results/{sample}/filter_host_reads/{sample}_dehuman.cram",
        checksum = "results/{sample}/filter_host_reads/{sample}_dehuman.cram.md5",
    params:
        sort_tmp=temp("results/{sample}/filter_host_reads/sort.tmp"),
        ref = config["resources"]["reference"],
    log:
        outfile="logs/{sample}/filter_host_reads/filter_host_reads.out.log",
        errfile="logs/{sample}/filter_host_reads/filter_host_reads.err.log",
    benchmark:
        "logs/benchmark/filter_host_reads/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        """
        samtools view -@ {threads} \
                    -h {input.bam} \
            | grep -v "_host"
            | samtools view -@ {threads} \
                    -bh -o {output.bam}

        samtools view -@ {threads} \
                    -h {input.bam} \
            | grep -v "_host" \
            | samtools sort -@ {threads} \
                    -T {params.sort_tmp} \
                    --reference {params.ref} \
                    -O cram \
                    -o {output.cram} \ 
                    2> >(tee -a {log.errfile} >&2)

        md5sum {output.cram} > {output.checksum} 2> >(tee -a {log.errfile} >&2)
        """