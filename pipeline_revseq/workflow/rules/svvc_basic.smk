rule fetch_primers:
    input:
        primerfile = config["resources"]["primer_file"],
        reference = config["resources"]["reference"]
    output:
        primerout = config["tools"]["fetch_primers"]["primerout"]
    benchmark:
        "logs/benchmark/fetch_primers/fetch_primers.benchmark"
    conda:
        "../envs/fetch_primers.yaml"
    shell:
        "python3 ../scripts/fetch_primers.py --input {input.primerfile} --reference {input.reference} --output {output.primerout}"


rule bwa_index:
    input:
        reference = config["resources"]["reference"]
    output:
        referenceout = config["resources"]["reference"] + '.bwt'
    benchmark:
        "logs/benchmark/bwa_index/bwa_index.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa index {input.reference}"


rule trim_galore:
    input:
        r1 = config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_1.fastq.gz",
        r2 = config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_2.fastq.gz",
    output:
        r1 = "results/{sample}/trim_galore/trimmed_r1.fq.gz",
        r2 = "results/{sample}/trim_galore/trimmed_r2.fq.gz",
    params:
        base_name = "data/{sample}",
        outdir = 'reuslts/{sample}/trim_galore',
        min_length = config["tools"]["trim_galore"]["min_length"],
        min_length_single = config["tools"]["trim_galore"]["min_length_single"],
    benchmark:
        "logs/benchmark/trim_galore/{sample}.benchmark"
    conda:
        "../envs/trim_galore.yaml"
    shell:
        'trim_galore --length {params.min_length} --output "{params.outdir}" --retain_unpaired --paired -r1 {params.min_length_single} -r2 {params.min_length_single} {input.r1} {input.r2} ; '
        "mv {params.outdir}/{sample}_1_val_1.fq.gz  {output.r1} ; "
        "mv {params.outdir}/{sample}_2_val_2.fq.gz  {output.r2}"


rule bwa:
    input:
        ref = config["resources"]["reference"],
        index = rules.bwa_index.output.referenceout,
        reads = rules.trim_galore.output
    output:
        outfile = "results/{sample}/bwa/mapped_reads.bam"
    benchmark:
        "logs/benchmark/bwa/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} |samtools view -Sb - > {output.outfile}"


rule fasta_index:
    input:
        ref = config["resources"]["reference"],
        primer_file = config["tools"]["fetch_primers"]["primerout"]
    output:
        outfile = "primers_regions.fasta",
        outindex = "primers_regions.fasta.faidx"
    benchmark:
        "logs/benchmark/fasta_index/fasta_index.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools faidx {input.ref} - o {output.outfile} --fai-idx {output.outindex}"


rule pileup:
    input:
        bam = rules.bwa.output.outfile,
        fasta = rules.fasta_index.output.outfile,
        primers = rules.fetch_primers.output.primerout
    output:
        outpile = "results/{sample}/pileup/allele_counts.npz"
    benchmark:
        "logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools mpileup -b {input.bam} -f {input.fasta} -l {input.primers} -o {output.outpile}"


rule pair_frequencies:
    input:
        alignment = rules.bwa.output.outfile
    output:
        outcounts = "results/{sample}/pair_freqs/pair_counts.pkl.gz"
    params:
        out_dir = "results/{sample}/pair_freqs"
    benchmark:
        "logs/benchmark/pair_frequencies/{sample}.benchmark"
    conda:
        "../envs/pair_frequencies.yaml"
    shell:
        "python3 ../scripts/pair_statistics.py --bam_file {input.alignment} --out_dir {params.out_dir}"


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
        min_cov = config["tools"]["consensus"]["min_cov"]
    benchmark:
        "logs/benchmark/consensus/{sample}.benchmark"
    conda:
        "../envs/consensus.yaml"
    shell:
        """
        python3 ../scripts/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} &
        python3 ../scripts/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} --min_freq {params.min_freq} --min_cov {params.min_cov}
        """

