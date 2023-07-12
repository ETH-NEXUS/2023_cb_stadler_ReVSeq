## Only this rule must be aware of the real sample names, nothing else
metadata_path = pd.read_table(config["resources"]["metadata_name"], header=None)
samples_table = pd.read_csv(metadata_path, sep=";", header=0)
samples = samples_table.set_index("sample", drop=False)
sample_ids = samples_table["sample"].tolist()
validate(samples, "../schema/sample_map.schema.yaml")
lanes = samples_table.set_index("lane", drop=False)
lane_ids = samples_table["lane"].tolist()

rule anonymize_samples:
    input:
        r1 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_{lane}_R1.fastq.gz", lane=lane_ids, sample = sample_ids),
        r2 = expand(config["inputOutput"]["input_fastqs"]+"/{{sample}}/{{sample}}_{lane}_R2.fastq.gz", lane=lane_ids, sample = sample_ids)
    output:
        outr1 = "results/{sample}/anonymize_samples/{ethid}_{lane}_R1.fastq.gz",
        outr2 = "results/{sample}/anonymize_samples/{ethid}_{lane}_R1.fastq.gz"
        outtable = "results/anonymize_samples/anonymization_table.tsv"
    params:
        sort_tmp=temp("results/{sample}/pileup/sort.tmp"),
    log:
        outfile="logs/{sample}/pileup/pileup.out.log",
        errfile="logs/{sample}/pileup/pileup.err.log",
    benchmark:
        "logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        python workflow/scripts/create_anonimization_table.py --r1 {input.r1} --r2 {input.r2}
        """

rule anonymize_samples:
    input:
        bam = rules.filter_host_reads.output.bam,
        fasta = rules.merge_refs.output.outref,
        #primers = config["resources"]["primer_file"],
    output:
        outpile = "results/{sample}/pileup/pileup.txt"
    params:
        sort_tmp=temp("results/{sample}/pileup/sort.tmp"),
    log:
        outfile="logs/{sample}/pileup/pileup.out.log",
        errfile="logs/{sample}/pileup/pileup.err.log",
    benchmark:
        "logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} \
                    -T {params.sort_tmp} \
                    -O bam {input.bam}
            | samtools mpileup -f {input.fasta} -o {output.outpile}
        """

rule validate_anonimization:
    input:
        bam = rules.filter_host_reads.output.bam,
        fasta = rules.merge_refs.output.outref,
        #primers = config["resources"]["primer_file"],
    output:
        outpile = "results/{sample}/pileup/pileup.txt"
    params:
        sort_tmp=temp("results/{sample}/pileup/sort.tmp"),
    log:
        outfile="logs/{sample}/pileup/pileup.out.log",
        errfile="logs/{sample}/pileup/pileup.err.log",
    benchmark:
        "logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} \
                    -T {params.sort_tmp} \
                    -O bam {input.bam}
            | samtools mpileup -f {input.fasta} -o {output.outpile}
        """