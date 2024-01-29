rule aggregate:
    input:
        assignment = expand(rules.assign_virus.output.table, sample=sample_ids),
        validation = expand(rules.validate_assignment.output.validation, sample=sample_ids),
        #bwa_all =expand(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_readcount_all.txt", sample=sample_ids),
        #bwa = expand(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/bwa/{sample}_readcount.txt",sample=sample_ids),
        #dehuman = expand(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_readcount.txt", sample=sample_ids),
        bwa_all = expand(rules.bwa.output.readcount_all, sample=sample_ids),
        bwa = expand(rules.bwa.output.readcount, sample=sample_ids),
        dehuman = expand(rules.dehuman.output.readcount, sample=sample_ids),
        filter_alignment = expand(rules.filter_alignment.output.readcount, sample=sample_ids),
        duplicate = expand(rules.remove_duplicates.output.readcount, sample=sample_ids),
    output:
        aggregated_assignment = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/aggregate/"+config["plate"]+"_aggregated_assignment_table.tsv",
        aggregated_qc = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/aggregate/"+config["plate"]+"_aggregated_qc_table.tsv",
    params:
        inputdir = config["inputOutput"]["output_dir"]+"/"+config["plate"],
        assigndirname = "assign_virus",
        validatedirname = "validate_assignment",
        bwaname = "bwa",
        dehumanname = "dehuman",
        filtername = "filter_alignment",
        duplicatename = "remove_duplicates",
        qcname = "qualimap_filtered",
        sample_map = config["sample_map"],
        pseudoanon_metadata = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_metadata.csv",
        match_table = config["tools"]["validate_assignment"]["match_table"],
        readnum_threshold = config['tools']['general']['min_readcount'],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/aggregate/aggregate.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/aggregate/aggregate.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/aggregate/aggregate.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/aggregate_tables.py \
        --inputdir {params.inputdir} \
        --assignment_subdir {params.assigndirname} \
        --validation_subdir {params.validatedirname} \
        --outfile {output.aggregated_assignment} \
        --sample_map {params.sample_map} \
        --pseudoanon_metadata {params.pseudoanon_metadata} \
        --match_table {params.match_table} \
        --readnum_threshold {params.readnum_threshold}

        python workflow/scripts/aggregate_qc.py \
        --inputdir {params.inputdir} \
        --bwa_subdir {params.bwaname} \
        --dehuman_subdir {params.dehumanname} \
        --filter_subdir {params.filtername} \
        --duplicates_subdir {params.duplicatename} \
        --outfile {output.aggregated_qc} \
        --sample_map {params.sample_map}

        """

