rule package_results:
    input:
        assignment = rules.assign_virus.output.table,
        rawr1 =rules.merge_lanes.output.r1,
        rawr2 = rules.merge_lanes.output.r2,
        bam = rules.remove_duplicates.output.bam,
        dehuman_cram = rules.cram.output.cram,
        fastqcr1 = rules.fastqc_merged.output.zip1,
        fastqcr2 = rules.fastqc_merged.output.zip2,
        multiqc = rules.multiqc.output.outfile,
    output:
        assignment = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/assign_virus/{sample}_count_table.tsv",
        rawr1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/raw_data/{sample}_merged_R1.fastq.gz",
        rawr2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/raw_data/{sample}_merged_R2.fastq.gz",
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/bam/{sample}_remove_duplicates.bam",
        dehuman_cram = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"package_results/cram/{sample}.cram",
        fastqcr1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"package_results/fastqc/{sample}_merged_R1.fastqc.zip",
        fastqcr2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"package_results/fastqc/{sample}_merged_R2.fastqc.zip",
    params:
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/",
        metadata = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        multiqcdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/multiqc",
        empty_samples = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        metadata_out = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"package_results/"+config["plate"]+"_metadata.csv",
        multiqc_out = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/multiqc",
        empty_samples_out = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"package_results/"+config["plate"]+"empty_samples.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/package_results/{sample}package_results.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/package_results/{sample}package_results.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/package_results/{sample}package_results.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        ln -s {input.assignment} {output.assignment}
        ln -s {input.rawr1} {output.rawr1}
        ln -s {input.rawr2} {output.rawr2}
        ln -s {input.bam} {output.bam}
        ln -s {input.dehuman_cram} {output.dehuman_cram}
        ln -s {input.fastqcr1} {output.fastqcr1}
        ln -s {input.fastqcr2} {output.fastqcr2}
		ln -s {params.multiqcdir} {params.multiqc_out}
		ln -s {params.metadata} {params.metadata_out}
		ln -s {params.empty_samples} {params.empty_samples_out}
        """