rule filter_alignment:
    input:
        bam = rules.dehuman.output.bam
    output:
        bam = (config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/filter_alignment/{sample}_filter_alignment.bam")#temp
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
        samtools view -bh -F 2304 -q 1 -e '![XA]' {input.bam} - | samtools collate -u -O {input.bam} - | samtools fixmate -u - | samtools view -f 1 -o out.bam
        """

rule sort:
    input:
        bwa = rules.filter_alignment.output.bam
    output:
        sorted = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/sort/{sample}_sorted_reads.bam",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/sort/sort.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/sort/sort.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/sort/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} --output-fmt=BAM -o {output.sorted} {input.bwa} 2> >(tee {log.errfile} >&2)"


rule remove_duplicates:
    input:
        bam = rules.sort.output.sorted
    output:
        bam = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/remove_duplicates/{sample}_remove_duplicates.bam",
        metrics = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/remove_duplicates/{sample}_metrics.txt"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/remove_duplicates/remove_duplicates.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/{sample}/remove_duplicates/remove_duplicates.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/logs/benchmark/remove_duplicates/{sample}.benchmark"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true 2> >(tee {log.errfile} >&2)"


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
        counts = rules.idxstats.output.counts
    output:
        table = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_count_table.tsv",
    params:
        ref_table = config["resources"]["reference_table"],
        prefix = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/assign_virus/{sample}_",
        outlier_percentile = config["tools"]["assign_virus"]["outlier_percentile"],
        outlier_percentile_collapsed = config["tools"]["assign_virus"]["outlier_percentile_collapsed"],
        lookup = config["tools"]["assign_virus"]["lookup"],
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
        --lookup {params.lookup}  2> >(tee {log.errfile} >&2)
        """


rule validate_assignment:
    input:
        assignment = rules.assign_virus.output.table,
    output:
        validation = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"/{sample}/validate_assignment/{sample}_validation.txt",
    params:
        metadata_dir = config["resources"]["metadata_dir"],
        pseudoanontable = config["inputOutput"]["output_dir"]+"/"+config["plate"]+"_pseudoanon_table.tsv",
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
        --anonymization_table {params.pseudoanontable} \
        --match_table {params.match_table} \
        --output {output.validation}  2> >(tee {log.errfile} >&2)
        """
