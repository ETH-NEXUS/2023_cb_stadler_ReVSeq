rule filter_alignment:
    input:
        #bam = rules.dehuman.output.bam
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/dehuman/{sample}_dehuman.bam"
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/filter_alignment/{sample}_filter_alignment.bam"),
        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/filter_alignment/{sample}_readcount.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/filter_alignment/filter_alignment.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/filter_alignment/filter_alignment.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/filter_alignment/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        # remove flags 256 and 2048, remove MAPQ = 0 and remove any read with tag XA
        # fix the mate pair flags to know which reads became unpaired after filtering
        # remove unpaired reads
        # https://github.com/samtools/samtools/issues/1668
        """
        samtools view -bh -F 2304 -q 1 -e '![XA]' {input.bam} | \
        samtools collate -u -O - | \
        samtools fixmate -u - - | \
        samtools view -h -f 1 | \
        samtools sort -@ {threads} --output-fmt=BAM -o {output.bam} - 2> >(tee {log.errfile} >&2)
        samtools view -c -F 4 {output.bam} > {output.readcount}
        """


rule remove_duplicates:
    input:
        bam = rules.filter_alignment.output.bam
    output:
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/remove_duplicates/{sample}_remove_duplicates.bam",
        metrics = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/remove_duplicates/{sample}_metrics.txt",
        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/remove_duplicates/{sample}_readcount.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/remove_duplicates/remove_duplicates.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/remove_duplicates/remove_duplicates.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/remove_duplicates/{sample}.benchmark"
    conda:
        "../envs/picard.yaml"
    shell:
        """
        picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true 2> >(tee {log.errfile} >&2)
        samtools view -c -F 4 {output.bam} > {output.readcount}
        """


rule qualimap_filtered:
    input:
        bam = rules.remove_duplicates.output.bam,
        readcount = rules.remove_duplicates.output.readcount,
    output:
        report = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap_filtered/qualimapReport.html",
        genome_res = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap_filtered/genome_results.txt",
        qc_status = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap_filtered/qc_status.txt",
    params:
        regions = config["resources"]["reference_table"],
        outdir = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/qualimap_filtered/",
        min_reads = config["tools"]["general"]["min_readcount"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/qualimap_filtered/qualimap.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/qualimap_filtered/qualimap.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/qualimap_filtered/qualimap.benchmark"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        if (( $(cat {input.readcount}) < {params.min_reads} ))
        then
            touch {output.report} {output.genome_res}
            echo "FAIL" | tee {output.qc_status}
        else
            echo "PASS" | tee {output.qc_status}
            unset DISPLAY
            qualimap bamqc -outdir {params.outdir} -bam {input.bam} --feature-file {params.regions} -c 2> >(tee {log.errfile} >&2)
        fi
        """


rule samtools_index:
    input:
        bam = rules.remove_duplicates.output.bam
    output:
        index = rules.remove_duplicates.output.bam + '.bai'
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/sort/samtools_index.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/sort/samtools_index.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/{sample}/sort/samtools_index.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bam} 2> >(tee {log.errfile} >&2)"



rule pileup:
    input:
        bam = rules.remove_duplicates.output.bam,
        fasta = config["resources"]["reference"]
    output:
        outpile = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/pileup/{sample}_pileup.txt"
    params:
        min_map_qual=config["tools"]["pileup"]["min_map_qual"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/pileup/pileup.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/pileup/pileup.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools mpileup -q {params.min_map_qual} -f {input.fasta} {input.bam} > {output.outpile}  2> >(tee {log.errfile} >&2)
        """


rule idxstats:
    input:
        bam = rules.remove_duplicates.output.bam,
        index = rules.samtools_index.output.index,
    output:
        idxstats = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/idxstats/{sample}_idxstats.txt",
        counts = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/idxstats/{sample}_counts.txt",
    params:
        min_map_qual=config["tools"]["idxstats"]["min_map_qual"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/idxstats/{sample}_idxstats.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/idxstats/{sample}_idxstats.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/idxstats/{sample}_idxstats.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -c -q {params.min_map_qual} -F 4 {input.bam} > {output.counts} 2> >(tee {log.errfile} >&2)
        samtools idxstats {input.bam} > {output.idxstats}  2> >(tee {log.errfile} >&2)
        """


rule assign_virus:
    input:
        idxstats = rules.idxstats.output.idxstats,
        counts = rules.idxstats.output.counts,
        genome_res = rules.qualimap_filtered.output.genome_res,
    output:
        table = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_substrain_count_table.tsv",
        strain_table = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_strain_count_table.tsv",
        boxplot = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_substrain_proportions_boxplot.pdf",
        strain_boxplot = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_strain_proportions_boxplot.pdf",
    params:
        ref_table = config["resources"]["reference_table"],
        prefix = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_",
        outlier_percentile = config["tools"]["assign_virus"]["outlier_percentile"],
        outlier_percentile_collapsed = config["tools"]["assign_virus"]["outlier_percentile_collapsed"],
        lookup = config["tools"]["assign_virus"]["lookup"],
        dp_threshold = config["tools"]["assign_virus"]["dp_threshold"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/assign_virus/{sample}_assignment.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/assign_virus/{sample}_assignment.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/assign_virus/{sample}_assignment.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
            python workflow/scripts/assign_virus.py \
            --ref_table {params.ref_table} \
            --idxstats {input.idxstats} \
            --counts {input.counts} \
            --out_prefix {params.prefix} \
            --outlier_percentile {params.outlier_percentile} \
            --outlier_percentile_collapsed {params.outlier_percentile_collapsed} \
            --lookup {params.lookup}  \
            --genome_res {input.genome_res} \
            --dp_threshold {params.dp_threshold} 2> >(tee {log.errfile} >&2)
        """


rule validate_assignment:
    input:
        assignment = rules.assign_virus.output.strain_table,
    output:
        validation = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/validate_assignment/{sample}_validation.tsv",
    params:
        metadata_dir = config["resources"]["metadata_dir"],
        pseudoanontable = config["inputOutput"]["input_fastqs"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv",
        match_table = config["tools"]["validate_assignment"]["match_table"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/validate_assignment/{sample}_validation.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/validate_assignment/{sample}_validation.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/validate_assignment/{sample}_validation.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
            python workflow/scripts/validate_assignment.py \
            --metadata_dir {params.metadata_dir} \
            --pseudoanon_table {params.pseudoanontable} \
            --ethid {wildcards.sample} \
            --match_table {params.match_table} \
            --count_table {input.assignment} \
            --output {output.validation}  2> >(tee {log.errfile} >&2)
        """


rule consensus:
    input:
        bam = rules.remove_duplicates.output.bam,
        qualimap_qc = rules.qualimap_filtered.output.qc_status,
    output:
        vcf = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_calls.vcf.gz",
        calls_norm = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_calls_norm.bcf",
        calls_norm_filt_indels = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_calls_norm_filt_indel.bcf",
        all_consensus = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_all_consensus.fa",
    params:
        ref = config["resources"]["reference"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/consensus/{sample}_consensus.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/consensus/{sample}_consensus.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/consensus/{sample}_consensus.benchmark"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        if [ $(cat {input.qualimap_qc}) = "PASS" ]
        then
            bcftools mpileup -Ou -f {params.ref} {input.bam} | bcftools call -mv -Oz -o {output.vcf}
            bcftools index {output.vcf}

            # normalize indels
            bcftools norm -f {params.ref} {output.vcf} -Ob -o {output.calls_norm}

            # filter adjacent indels within 5bp
            bcftools filter --IndelGap 5 {output.calls_norm} -Ob -o {output.calls_norm_filt_indels}
            bcftools index {output.calls_norm_filt_indels}

            # apply variants to create consensus sequence
            cat {params.ref} | bcftools consensus {output.calls_norm_filt_indels} > {output.all_consensus}
        else
            touch {output.vcf} {output.calls_norm} {output.calls_norm_filt_indels} {output.all_consensus}
        fi
        """


rule postprocess_consensus:
    input:
        all_consensus = rules.consensus.output.all_consensus,
        assignment = rules.assign_virus.output.strain_table,
        qualimap_qc = rules.qualimap_filtered.output.qc_status,
    output:
        consensus = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/postprocess_consensus/{sample}_consensus.fa",
    params:
        ref = config["resources"]["reference_table"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/postprocess_consensus/{sample}_consensus.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/postprocess_consensus/{sample}_consensus.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/postprocess_consensus/{sample}_consensus.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        if [ $(cat {input.qualimap_qc}) = "PASS" ]
        then
            # Fetch only the consensus in the regions of the most represented virus
            python workflow/scripts/filter_consensus.py \
            --ref_table {params.ref} \
            --assignment {input.assignment} \
            --consensus {input.all_consensus} \
            --output {output.consensus}  2> >(tee {log.errfile} >&2)
        else
            touch {output.consensus}
        fi
        """
