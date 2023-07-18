rule pileup:
    input:
        bam = rules.filter_host_reads.output.bam,
        fasta = config["resources"]["reference"]
        #fasta = rules.merge_refs.output.outref,
        #primers = config["resources"]["primer_file"],
    output:
        outpile = config["inputOutput"]["output_dir"]+"/{sample}/pileup/pileup.txt"
    params:
        sort_tmp=(config["inputOutput"]["output_dir"]+"/{sample}/pileup/sort.tmp"),
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/pileup/pileup.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/pileup/pileup.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools sort -@ {threads} \
                    -T {params.sort_tmp} \
                    -O bam {input.bam}
            | samtools mpileup -f {input.fasta} -o {output.outpile}
        """


rule pair_frequencies:
    input:
        alignment = rules.bwa.output.bam
    output:
        outcounts = config["inputOutput"]["output_dir"]+"/{sample}/pair_freqs/pair_counts.pkl.gz"
    params:
        out_dir = config["inputOutput"]["output_dir"]+"/{sample}/pair_freqs"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/pair_freqs/pair_freqs.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/pair_freqs/pair_freqs.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/pair_frequencies/{sample}.benchmark"
    conda:
        "../envs/pair_frequencies.yaml"
    shell:
        "python workflow/scripts/pair_statistics.py --bam_file {input.alignment} --out_dir {params.out_dir}"


rule consensus:
    input:
        counts = rules.pileup.output.outpile
    output:
        config["inputOutput"]["output_dir"]+"/{sample}/consensus/consensus.fasta",
        config["inputOutput"]["output_dir"]+"/{sample}/consensus/figures/coverage.png",
        config["inputOutput"]["output_dir"]+"/{sample}/consensus/figures/diversity.png",
        config["inputOutput"]["output_dir"]+"/{sample}/consensus/minor.fasta"
    params:
        out_dir = config["inputOutput"]["output_dir"]+"/{sample}/consensus",
        min_freq = config["tools"]["consensus"]["min_freq"],
        min_cov = config["tools"]["consensus"]["min_cov"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/consensus/consensus.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/consensus/consensus.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/consensus/{sample}.benchmark"
    conda:
        "../envs/consensus.yaml"
    shell:
        """
        python workflow/scripts/coverage_consensus_diversity.py --sample {params.out_dir} --out_dir {params.out_dir} &
        python workflow/scripts/minor_variant.py --sample {params.out_dir} --out_dir {params.out_dir} --min_freq {params.min_freq} --min_cov {params.min_cov}
        """

