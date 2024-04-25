rule filter_alignment:
    input:
        bam = rules.dehuman.output.bam
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/filter_alignment/{sample}_filter_alignment.bam"),
        readcount = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/filter_alignment/{sample}_readcount.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/filter_alignment/filter_alignment.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/filter_alignment/filter_alignment.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/filter_alignment/{sample}_filter_alignment.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
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
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/remove_duplicates/{sample}_remove_duplicate.benchmark"
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
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/qualimap_filtered/{sample}_qualimap_filtered.benchmark"
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
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/sort/{sample}_samtools_index.benchmark"
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
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/pileup/{sample}_pileup.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools mpileup -q {params.min_map_qual} -f {input.fasta} {input.bam} > {output.outpile}  2> >(tee {log.errfile} >&2)
        """

rule depth:
    input:
        bam = rules.remove_duplicates.output.bam,
        index = rules.samtools_index.output.index,
    output:
        depth = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/depth/{sample}_depth.tsv",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/depth/{sample}_depth.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/depth/{sample}_depth.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/depth/{sample}_depth.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools depth -a -o {output.depth} {input.bam}
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
        depth = rules.depth.output.depth,
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
        lookup = config["tools"]["general"]["lookup"],
        coverage_threshold = config["tools"]["assign_virus"]["coverage_threshold"],
        readnum_threshold = config["tools"]["general"]["min_readcount"],
        taxon = config["tools"]["assign_virus"]["taxon"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/assign_virus/{sample}_assignment.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/assign_virus/{sample}_assignment.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/assign_virus/{sample}_assign_virus.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
            python workflow/scripts/assign_virus.py \
            --ref_table {params.ref_table} \
            --idxstats {input.idxstats} \
            --counts {input.counts} \
            --taxon_table {params.taxon} \
            --depth_table {input.depth} \
            --out_prefix {params.prefix} \
            --outlier_percentile {params.outlier_percentile} \
            --outlier_percentile_collapsed {params.outlier_percentile_collapsed} \
            --lookup {params.lookup}  \
            --genome_res {input.genome_res} \
            --coverage_threshold {params.coverage_threshold} \
            --readnum_threshold {params.readnum_threshold} 2> >(tee {log.errfile} >&2)
        """


rule validate_assignment:
    input:
        assignment = rules.assign_virus.output.strain_table,
    output:
        validation = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/validate_assignment/{sample}_validation.tsv",
    params:
        metadata_dir = config["resources"]["metadata_dir"],
        pseudoanontable = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/"+config["plate"]+"_pseudoanon_table.tsv",
        match_table = config["tools"]["validate_assignment"]["match_table"],
        ctrl_pos = config["tools"]["general"]["pos"],
        ctrl_neg = config["tools"]["general"]["neg"],
        plate = config["plate"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/validate_assignment/{sample}_validation.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/validate_assignment/{sample}_validation.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/validate_assignment/{sample}_validate_assignment.benchmark"
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
            --controls {params.ctrl_pos},{params.ctrl_neg} \
            --plate {params.plate} \
            --output {output.validation}  2> >(tee {log.errfile} >&2)
        """

#rule create_chromosome_file:
#    input:
#        assignment = rules.assign_virus.output.table,
#    output:
#        chrom_file = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/create_chromosome_file/{sample}_chromosome_file.txt",
#        chrom_file_gzip = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/create_chromosome_file/{sample}_chromosome_file.txt.gz",
#    log:
#        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/create_chromosome_file/{sample}_create_chromosome_file.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/create_chromosome_file/{sample}_create_chromosome_file.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/create_chromosome_file/{sample}_create_chromosome_file.benchmark"
#    shell:
#        """
#        name=sort -r -t$'\t' -k6 {input.assignment} | head -n 2 | tail -n 1 | awk '{print $1}'
#        echo "${name}\t1\t"
#        """"

rule consensus:
    input:
        bam = rules.remove_duplicates.output.bam,
        #qualimap_qc = rules.qualimap_filtered.output.qc_status,
    output:
        calls_norm_filt_indels = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_calls_norm_filt_indel.bcf",
        all_consensus = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_all_consensus.fa",
        low_cov_bed = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/consensus/{sample}_low_coverage_sites.bed",
    params:
        ref = config["resources"]["reference"],
        mincov = config["tools"]["consensus"]["minimum_coverage"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/consensus/{sample}_consensus.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/consensus/{sample}_consensus.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/consensus/{sample}_consensus.benchmark"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
            bcftools mpileup -Ou -f {params.ref} {input.bam} | \
            bcftools call --ploidy 1 -mv -Ou | \
            bcftools norm -f {params.ref} -Ou | \
            bcftools filter --IndelGap 5 -Ob -o {output.calls_norm_filt_indels}
            bcftools index -o {output.calls_norm_filt_indels}.csi {output.calls_norm_filt_indels}

            bedtools genomecov -ibam {input.bam} -bga | awk '$4 < {params.mincov}' > {output.low_cov_bed}

            # apply variants to create consensus sequence
            bcftools consensus {output.calls_norm_filt_indels} -f {params.ref} -m {output.low_cov_bed} > {output.all_consensus}
        """


rule postprocess_consensus:
    input:
        all_consensus = rules.consensus.output.all_consensus,
        assignment = rules.assign_virus.output.table,
        qualimap_qc = rules.qualimap_filtered.output.qc_status,
    output:
        consensus = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/postprocess_consensus/{sample}_consensus.fa",
        consensus_gzip = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/postprocess_consensus/{sample}_consensus.fa.gz",
        count_n = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/postprocess_consensus/{sample}_count_n.txt",
    params:
        ref = config["resources"]["reference_table"],
        consensus_type = config["tools"]["consensus"]["consensus_type"],
        lookup = config["tools"]["general"]["lookup"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/postprocess_consensus/{sample}_consensus.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/postprocess_consensus/{sample}_consensus.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/postprocess_consensus/{sample}_postprocess_consensus.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        # Fetch only the consensus in the regions of the most represented virus
        python workflow/scripts/filter_consensus.py \
        --ref_table {params.ref} \
        --assignment {input.assignment} \
        --consensus {input.all_consensus} \
        --consensus_type {params.consensus_type} \
        --output {output.consensus} 2> >(tee {log.errfile} >&2)
        
        gzip -c {output.consensus} > {output.consensus_gzip}

        python workflow/scripts/count_n.py \
        --input {output.consensus} \
        --output {output.count_n}
        """
