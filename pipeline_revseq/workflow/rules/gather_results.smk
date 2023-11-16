rule gather_results_samples:
    input:
        assignment = rules.assign_virus.output.table,
        rawr1 = rules.merge_lanes.output.r1,
        rawr2 = rules.merge_lanes.output.r2,
        bam = rules.remove_duplicates.output.bam,
        dehuman_cram = rules.cram.output.cram,
        fastqcr1 = rules.fastqc_merged.output.zip1,
        fastqcr2 = rules.fastqc_merged.output.zip2,
        multiqc = rules.multiqc.output.outfile,
        consensus = rules.postprocess_consensus.output.consensus,
        qc_status = rules.qualimap_filtered.output.qc_status,
    output:
        ####to_upload = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/{sample}/{sample}_to_upload.txt",
        assignment = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_count_table.tsv",
        rawr1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_merged_R1.fastq.gz",
        rawr2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_merged_R2.fastq.gz",
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_remove_duplicates.bam",
        dehuman_cram = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}.cram",
        fastqcr1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_merged_R1.fastqc.zip",
        fastqcr2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_R2.fastqc.zip",
        consensus = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/sample_{sample}/{sample}_consensus.fa",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/gather_results_sample/{sample}_gather_results.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/gather_results_sample/{sample}_gather_results.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/gather_results_sample/{sample}_gather_results.benchmark"
    shell:
        """
        ln -s {input.assignment} {output.assignment}
        ln -s {input.rawr1} {output.rawr1}
        ln -s {input.rawr2} {output.rawr2}
        ln -s {input.bam} {output.bam}
        ln -s {input.dehuman_cram} {output.dehuman_cram}
        ln -s {input.fastqcr1} {output.fastqcr1}
        ln -s {input.fastqcr2} {output.fastqcr2}
        ln -s {input.consensus} {output.consensus}
        """

rule gather_results_plate:
    input:
        metadata = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        multiqcdir = rules.multiqc.output.outdir,
        multiqcdir_filtered = rules.multiqc_filtered.output.outdir,
        empty_samples = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        version = config["resources"]["pipeline_version"],
    output:
        metadata = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        multiqcdir = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/multiqc"),
        multiqcdir_filtered = directory(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/multiqc_filtered"),
        empty_samples = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        version = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/"+config["plate"]+"/"+config["plate"]+"_pipeline_version.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/gather_results_plate/gather_results_plate.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/gather_results_plate/gather_results_plate.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/gather_results/gather_results_plate.benchmark"
    shell:
        """
        ln -s {input.metadata} {output.metadata}
        mkdir {output.multiqcdir}
        ln -s {input.multiqcdir}/multiqc_report.html {output.multiqcdir}
        mkdir {output.multiqcdir}/multiqc_data
        ln -s {input.multiqcdir}/multiqc_data/* {output.multiqcdir}/multiqc_data
        mkdir {output.multiqcdir_filtered}
        ln -s {input.multiqcdir_filtered}/multiqc_report.html {output.multiqcdir_filtered}
        mkdir {output.multiqcdir_filtered}/multiqc_data
        ln -s {input.multiqcdir_filtered}/multiqc_data/* {output.multiqcdir_filtered}/multiqc_data
        ln -s {input.empty_samples} {output.empty_samples}
        cat {input.version} > {output.version}
        """

#rule push_to_db:
#    input:
#        file_list = rules.gather_results,
#    output:
#        db_upload_success = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/push_to_db/db_upload_success",
#    log:
#		outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/push_to_db/push_to_db.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/push_to_db/push_to_db.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/push_to_db/push_to_db.benchmark"
#    shell:
#        """
#        (XXXX && \
#        echo "SUCCESS" > {output.db_upload_response}) || \
#        exit 1
#		"""
