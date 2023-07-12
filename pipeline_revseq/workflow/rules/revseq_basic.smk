#rule sort:
#    input:
#        bwa = rules.filter_host_reads.output.outbam
#    output:
#        outfile = tmp("results/{sample}/filter_host_reads/sorted_reads.bam")
#    log:
#        outfile="logs/{sample}/sort/sort.out.log",
#        errfile="logs/{sample}/sort/sort.err.log",
#    benchmark:
#        "logs/benchmark/sort/{sample}.benchmark"
#    conda:
#        "../envs/bwa.yaml"
#    threads: config["threads"]
#    shell:
#        "samtools sort -@ {threads} --output-fmt=BAM -o {output.outfile} {input.bwa}"


rule pileup:
    input:
        bam = rules.filter_host_reads.output.bam,
        fasta = config["resources"]["reference"]
        #fasta = rules.merge_refs.output.outref,
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


rule pair_frequencies:
    input:
        alignment = rules.bwa.output.bam
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

