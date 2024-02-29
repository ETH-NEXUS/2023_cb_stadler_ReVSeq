rule gather_results_samples:
    input:
        assignment = rules.assign_virus.output.table,
        rawr1 = rules.merge_lanes.output.r1,
        rawr2 = rules.merge_lanes.output.r2,
        bam = rules.remove_duplicates.output.bam,
        dehuman_cram = rules.cram.output.cram,
        fastqcr1 = rules.fastqc_merged.output.zip1,
        fastqcr2 = rules.fastqc_merged.output.zip2,
        consensus = rules.postprocess_consensus.output.consensus,
        qc_status = rules.qualimap_filtered.output.qc_status,
    output:
        assignment = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_count_table.tsv",
        rawr1 = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_merged_R1.fastq.gz",
        rawr2 = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_merged_R2.fastq.gz",
        bam = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_remove_duplicates.bam",
        dehuman_cram = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}.cram",
        fastqcr1 = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_merged_R1.fastqc.zip",
        fastqcr2 = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_R2.fastqc.zip",
        consensus = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_consensus.fa",
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
        aggregated_assignment = rules.aggregate.output.aggregated_assignment,
        aggregated_qc = rules.aggregate.output.aggregated_qc,
    output:
        metadata = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        multiqcdir = directory(config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/multiqc"),
        multiqcdir_filtered = directory(config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/multiqc_filtered"),
        empty_samples = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        version = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_pipeline_version.txt",
        aggregated_assignment = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_aggregated_assignment.tsv",
        aggregated_qc = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_aggregated_qc.tsv",
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
        ln -s {input.aggregated_assignment} {output.aggregated_assignment}
        ln -s {input.aggregated_qc} {output.aggregated_qc}
        cat {input.version} > {output.version}
        """

envvars:
    "USERNAME_REVSEQ",
    "PASSWORD_REVSEQ",
rule push_to_db:
    input:
        multiqcdir_filtered = rules.gather_results_plate.output.multiqcdir_filtered,
    output:
        db_upload_status = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/db_upload_status",
    params:
        #credentials_file = config["tools"]["push_to_db"]["credentials_file"],
        plate = config["plate"],
        revseq_username = os.environ["USERNAME_REVSEQ"],
        revseq_password = os.environ["PASSWORD_REVSEQ"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/push_to_db/push_to_db.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/push_to_db/push_to_db.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/push_to_db/push_to_db.benchmark"
    conda:
        "../envs/revseqdataloader.yaml"
    shell:
        """
        echo "Starting db_upload of plate {params.plate}"
        export USERNAME_REVSEQ={params.revseq_username}
        export PASSWORD_REVSEQ={params.revseq_password}
        (python workflow/scripts/upload_to_db.py --plate {params.plate} && \
        touch {output.db_upload_status}) || \
        echo "Failed uploading to database"
		"""

rule viollier_upload:
    input:
        gathered_results = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/gather_results/",
    output:
        viollier_upload_success = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/viollier_upload/viollier_upload_success",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/viollier_upload/viollier_upload.benchmark",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/viollier_upload/viollier_upload.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/viollier_upload/viollier_upload.err.log",
    params:
        revseq_executable = config["tools"]["viollier_upload"]["revseq_executable"],
    shell:
        """
        ( {params.revseq_executable} uploadviollier {input.gathered_results} && \
        echo "SUCCESS" > {output.viollier_upload_success}) || \
        exit 1
		"""