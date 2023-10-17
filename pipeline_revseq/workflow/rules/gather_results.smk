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
        to_upload = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_to_upload.txt",
        #assignment = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_count_table.tsv",
        #rawr1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_merged_R1.fastq.gz",
        #rawr2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_merged_R2.fastq.gz",
        #bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_remove_duplicates.bam",
        #dehuman_cram = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}.cram",
        #fastqcr1 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_merged_R1.fastqc.zip",
        #fastqcr2 = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_R2.fastqc.zip",
        #consensus = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/{sample}/{sample}_consensus.fa",
    params:
        metadata = config["inputOutput"]["input_fastqs"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        multiqcdir = directory(rules.multiqc.output.outdir),
        empty_samples = config["inputOutput"]["input_fastqs"]+"/"+config["plate"]+"/"+config["plate"]+"_empty_samples.txt",
        version = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/package_results/"+config["plate"]+"/"+config["plate"]+"_pipeline_version.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/package_results/{sample}_package_results.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/package_results/{sample}_package_results.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/package_results/{sample}_package_results.benchmark"
    shell:
        """
        echo {input.assignment} > {output.to_upload}
        echo {input.rawr1} >> {output.to_upload}
        echo {input.rawr2} >> {output.to_upload}
        echo {input.bam} >> {output.to_upload}
        echo {input.dehuman_cram} >> {output.to_upload}
        echo {input.fastqcr1} >> {output.to_upload}
        echo {input.fastqcr2} >> {output.to_upload}
        echo {input.consensus} >> {output.to_upload}
        echo {params.metadata} >> {output.to_upload}
        echo {params.empty_samples} >> {output.to_upload}
        echo {params.multiqcdir} >> {output.to_upload}
        echo {params.version} >> {output.to_upload}
        echo {input.qc_status} >> {output.to_upload}
        """


rule push_to_db:
    input:
        file_list = rules.gather_results,
    output:
        db_upload_response = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/push_to_db/db_upload_response.txt",
    log:
		outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/push_to_db/push_to_db.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/push_to_db/push_to_db.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/push_to_db/push_to_db.benchmark"
    shell:
        """
        (XXXX && \
        echo "SUCCESS" > {output.db_upload_response}) || \
        exit 1
		"""

rule upload_viollier:
    input:
        db_upload_response = rules.push_to_db.output.db_upload_response,
    output:
        viollier_upload_response = onfig["inputOutput"]["output_dir"]+"/"+config["plate"]+"/upload_viollier/viollier_upload_response.txt",
    log:
		outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/upload_viollier/upload_viollier.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/upload_viollier/upload_viollier.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/upload_viollier/upload_viollier.benchmark"
    conda:
        "../envs/viollier_upload.yaml"
    shell:
        """
        lftp XXXXX