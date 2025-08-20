rule gather_results_samples:
    input:
        assignment = rules.assign_virus.output.table,
        rawr1 = rules.merge_lanes.output.r1,
        rawr2 = rules.merge_lanes.output.r2,
        bam = rules.remove_duplicates.output.bam,
        dehuman_cram = rules.cram.output.cram,
        consensus = rules.postprocess_consensus.output.consensus,
        count_n = rules.postprocess_consensus.output.count_n,
        #consensus_cds = rules.consensus_cds.output.consensus_cds,
        #count_n_cds = rules.count_n_cds.output.count_n_cds,
        qc_status = rules.qualimap_filtered.output.qc_status,
        chr_file = rules.chr_file.output.chr_file,
        chr_file_gzip = rules.chr_file.output.chr_file_gzip,
        consensus_upload_gzip = rules.postprocess_consensus.output.consensus_upload_gzip,
    output:
        assignment = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_count_table.tsv",
        rawr1 = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_merged_R1.fastq.gz",
        rawr2 = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_merged_R2.fastq.gz",
        bam = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_remove_duplicates.bam",
        dehuman_cram = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}.cram",
        consensus = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_consensus.fa",
        complete = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/complete.txt",
        count_n = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_count_n.txt",
        #consensus_cds = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_consensus_cds.fa",
        #count_n_cds = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_count_n_cds.txt",
        chr_file = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/chr_file.txt",
        chr_file_gzip = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/chr_file.txt.gz",
        consensus_upload_gzip = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/sample_{sample}/{sample}_consensus_upload.gz",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/gather_results_sample/{sample}_gather_results.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/gather_results_sample/{sample}_gather_results.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/gather_results_sample/{sample}_gather_results.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
            python workflow/scripts/add_consensus_info.py --assignment {input.assignment} --consensus_n {input.count_n} --output {output.assignment}
            cp {input.rawr1} {output.rawr1}
            cp {input.rawr2} {output.rawr2}
            cp {input.bam} {output.bam}
            cp {input.consensus} {output.consensus}
            cp {input.count_n} {output.count_n}
            cp {input.dehuman_cram} {output.dehuman_cram}
            cp {input.chr_file} {output.chr_file}
            cp {input.chr_file_gzip} {output.chr_file_gzip}
            cp {input.consensus_upload_gzip} {output.consensus_upload_gzip}
            touch {output.complete}
        """

rule gather_results_plate:
    input:
        metadata = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        #multiqcdir = rules.multiqc.output.outdir,
        #multiqcdir_filtered = rules.multiqc_filtered.output.outdir,
        #multiqcdir_trimmed = rules.multiqc_trimmed.output.outdir,
        empty_samples = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        version = config["resources"]["pipeline_version"],
        #aggregated_assignment = rules.aggregate.output.aggregated_assignment,
        #aggregated_qc = rules.aggregate.output.aggregated_qc,
    output:
        metadata = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        #multiqcdir = directory(config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/multiqc"),
        #multiqcdir_filtered = directory(config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/multiqc_filtered"),
        #multiqcdir_trimmed = directory(config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/multiqc_trimmed"),
        empty_samples = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        version = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_pipeline_version.txt",
        #aggregated_assignment = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_aggregated_assignment.tsv",
        #aggregated_qc = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/"+config["plate"]+"_aggregated_qc.tsv",
        complete = config["tools"]["gather_results"]["outdir"]+"/"+config["plate"]+"/complete.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/gather_results_plate/gather_results_plate.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/gather_results_plate/gather_results_plate.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/gather_results/gather_results_plate.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/clean_metadata.py --input {input.metadata} --output {output.metadata}
        cp {input.empty_samples} {output.empty_samples}
        cat {input.version} > {output.version}
        touch {output.complete}
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
        samples_complete = expand(rules.gather_results_samples.output.complete, sample=sample_ids),
        plate_complete = expand(rules.gather_results_plate.output.complete, sample=sample_ids),
    output:
        viollier_upload_success = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/viollier_upload/viollier_upload_success",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/viollier_upload/viollier_upload.benchmark",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/viollier_upload/viollier_upload.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/viollier_upload/viollier_upload.err.log",
    params:
        revseq_executable = config["tools"]["viollier_upload"]["revseq_executable"],
        plate = config["plate"]
    shell:
        """
        ( {params.revseq_executable} uploadviollier {params.plate} && \
        echo "SUCCESS" > {output.viollier_upload_success}) || \
        exit 1
		"""