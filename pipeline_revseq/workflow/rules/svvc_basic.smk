rule fetch_primers:
    input:
        primerfile = config["resources"]["primer_file"],
        reference = config["resources"]["reference"]
    output:
        primerout = config["tools"]["fetch_primers"]["primerout"]
    log:
        outfile="logs/fetch_primers/fetch_primers.out.log",
        errfile="logs/fetch_primers/fetch_primers.err.log",
    benchmark:
        "logs/benchmark/fetch_primers/fetch_primers.benchmark"
    conda:
        "../envs/fetch_primers.yaml"
    shell:
        "bedtools getfasta -name -tab -fi {input.reference} -bed {input.primerfile} -fo {output.primerout}"
        #"python workflow/scripts/fetch_primers.py --primerfile {input.primerfile} --reference {input.reference} --output {output.primerout}"


rule bwa_index:
    input:
        reference = config["resources"]["reference"]
    output:
        referenceout = config["resources"]["reference"] + '.bwt'
    log:
        outfile="logs/bwa_index/bwa_index.out.log",
        errfile="logs/bwa_index/bwa_index.err.log",
    benchmark:
        "logs/benchmark/bwa_index/bwa_index.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.reference}"


rule trim_galore:
    input:
        r1 = config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_R1.fastq.gz",
        r2 = config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_R2.fastq.gz",
    output:
        r1 = "results/{sample}/trim_galore/trimmed_r1.fq.gz",
        r2 = "results/{sample}/trim_galore/trimmed_r2.fq.gz",
    params:
        base_name = "data/{sample}",
        outdir = 'results/{sample}/trim_galore',
        min_length = config["tools"]["trim_galore"]["min_length"],
        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
    log:
        outfile="logs/{sample}/trim_galore/trim_galore.out.log",
        errfile="logs/{sample}/trim_galore/trim_galore.err.log",
    benchmark:
        "logs/benchmark/trim_galore/{sample}.benchmark"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        'trim_galore --length {params.min_length} --output "{params.outdir}" --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {input.r1} {input.r2} ; '
        "mv {params.outdir}/{wildcards.sample}_R1_val_1.fq.gz  {output.r1} ; "
        "mv {params.outdir}/{wildcards.sample}_R2_val_2.fq.gz  {output.r2}"


rule bwa:
    input:
        ref = config["resources"]["reference"],
        index = rules.bwa_index.output.referenceout,
        reads = rules.trim_galore.output
    output:
        outfile = "results/{sample}/bwa/mapped_reads.bam"
    log:
        outfile="logs/{sample}/bwa/bwa.out.log",
        errfile="logs/{sample}/bwa/bwa.err.log",
    benchmark:
        "logs/benchmark/bwa/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} |samtools view -Sb - > {output.outfile}"


rule sort:
    input:
        bwa = rules.bwa.output.outfile
    output:
        outfile = "results/{sample}/sort/sorted_reads.bam"
    log:
        outfile="logs/{sample}/sort/sort.out.log",
        errfile="logs/{sample}/sort/sort.err.log",
    benchmark:
        "logs/benchmark/sort/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} --output-fmt=BAM -o {output.outfile} {input.bwa}"



rule fasta_index:
    input:
        ref = config["resources"]["reference"],
    output:
        outindex = config["resources"]["reference_dir"] + "primers_regions.fasta.faidx"
    log:
        outfile="logs/fasta_index/fasta_index.out.log",
        errfile="logs/fasta_index/fasta_index.err.log",
    benchmark:
        "logs/benchmark/fasta_index/fasta_index.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.ref} --fai-idx {output.outindex}"


rule pileup:
    input:
        bam = rules.sort.output.outfile,
        fasta = config["resources"]["reference"],
        primers = config["resources"]["primer_file"],
    output:
        outpile = "results/{sample}/pileup/allele_counts.npz"
    log:
        outfile="logs/{sample}/pileup/pileup.out.log",
        errfile="logs/{sample}/pileup/pileup.err.log",
    benchmark:
        "logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools mpileup -f {input.fasta} -l {input.primers} -o {output.outpile} {input.bam}"


rule pair_frequencies:
    input:
        alignment = rules.bwa.output.outfile
    output:
        outcounts = "results/{sample}/pair_freqs/pair_counts.pkl.gz"
    params:
        out_dir = "results/{sample}/pair_freqs"
    log:
        outfile="logs/{sample}/pair_freqs/pair_freqs.out.log",
        errfile="logs/{sample}/pair_freqs/pair_freqs.err.log",
    benchmark:
        "logs/benchmark/pair_frequencies/{sample}.benchmark"
    conda:
        "../envs/pair_frequencies.yaml"
    shell:
        "python workflow/scripts/pair_statistics.py --bam_file {input.alignment} --out_dir {params.out_dir}"


rule consensus:
    input:
        counts = rules.pileup.output.outpile
    output:
        "results/{sample}/consensus/consensus.fasta",
        "results/{sample}/consensus/figures/coverage.png",
        "results/{sample}/consensus/figures/diversity.png",
        "results/{sample}/consensus/minor.fasta"
    params:
        out_dir = "results/{sample}/consensus",
        min_freq = config["tools"]["consensus"]["min_freq"],
        min_cov = config["tools"]["consensus"]["min_cov"],
    log:
        outfile="logs/{sample}/consensus/consensus.out.log",
        errfile="logs/{sample}/consensus/consensus.err.log",
    benchmark:
        "logs/benchmark/consensus/{sample}.benchmark"
    conda:
        "../envs/consensus.yaml"
    shell:
        """
        python workflow/scripts/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} &
        python workflow/scripts/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} --min_freq {params.min_freq} --min_cov {params.min_cov}
        """

