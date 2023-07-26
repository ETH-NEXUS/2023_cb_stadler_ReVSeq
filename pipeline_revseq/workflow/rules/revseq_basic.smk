rule bwa:
    input:
        ref = config["resources"]["reference"],
        #ref = rules.merge_refs.output.outref,
        index = rules.bwa_index.output.referenceout,
        reads = rules.trim_galore.output
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/{sample}/bwa/{sample}_mapped_reads.bam")
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/bwa/bwa.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/bwa/bwa.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/bwa/{sample}.benchmark"
    conda:
        "../envs/bwa.yaml"
    shell:
        "bwa mem {input.ref} {input.reads} |samtools view -Sb - > {output.bam} 2> >(tee {log.errfile} >&2)"


rule remove_multimappers:
    input:
        bam = rules.bwa.output.bam
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/{sample}/remove_multimappers/{sample}_multimappers_removed.bam")
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/remove_multimappers/remove_multimappers.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/remove_multimappers/remove_multimappers.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/remove_multimappers/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools view -bh -F 256 {input.bam} -o {output.bam} 2> >(tee {log.errfile} >&2)"


rule remove_duplicates:
    input:
        bam = rules.remove_multimappers.output.bam
    output:
        bam = temp(config["inputOutput"]["output_dir"]+"/{sample}/remove_duplicates/{sample}_remove_duplicates.bam"),
        metrics = config["inputOutput"]["output_dir"]+"/{sample}/remove_duplicates/{sample}_metrics.txt"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/remove_duplicates/remove_duplicates.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/remove_duplicates/remove_duplicates.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/remove_duplicates/{sample}.benchmark"
    conda:
        "../envs/picard.yaml"
    shell:
        "picard MarkDuplicates I={input.bam} O={output.bam} M={output.metrics} REMOVE_DUPLICATES=true 2> >(tee {log.errfile} >&2)"


rule sort:
    input:
        bwa = rules.remove_duplicates.output.bam
    output:
        sorted = config["inputOutput"]["output_dir"]+"/{sample}/sort/{sample}_sorted_reads.bam",
        counts = config["inputOutput"]["output_dir"]+"/{sample}/sort/{sample}_counts.txt"
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/sort/sort.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/sort/sort.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/sort/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml"
    threads: config["threads"]
    shell:
        "samtools sort -@ {threads} --output-fmt=BAM -o {output.sorted} {input.bwa} 2> >(tee {log.errfile} >&2)"


rule samtools_index:
    input:
        bam = rules.sort.output.sorted
    output:
        index = rules.sort.output.sorted + '.bai'
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/sort/samtools_index.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/sort/samtools_index.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/{sample}/sort/samtools_index.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index {input.bam} 2> >(tee {log.errfile} >&2)"


rule pileup:
    input:
        bam = rules.sort.output.sorted,
        fasta = config["resources"]["reference"]
        #fasta = rules.merge_refs.output.outref,
        #primers = config["resources"]["primer_file"],
    output:
        outpile = config["inputOutput"]["output_dir"]+"/{sample}/pileup/{sample}_pileup.txt"
    params:
        sort_tmp=temp(config["inputOutput"]["output_dir"]+"/{sample}/pileup/sort.tmp"),
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/pileup/pileup.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/pileup/pileup.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/pileup/{sample}.benchmark"
    conda:
        "../envs/samtools.yaml 2> >(tee {log.errfile} >&2)"
    shell:
        """
        samtools sort -@ {threads} \
                    -T {params.sort_tmp} \
                    -O bam --verbosity 5 {input.bam} \
                | samtools mpileup -f {input.fasta} - > {output.outpile}  2> >(tee {log.errfile} >&2)
        """


rule idxstats:
    input:
        bam = rules.sort.output.sorted,
        index = rules.samtools_index.output.index
    output:
        idxstats = config["inputOutput"]["output_dir"]+"/{sample}/idxstats/{sample}_idxstats.txt",
        counts = config["inputOutput"]["output_dir"]+"/{sample}/idxstats/{sample}_counts.txt",
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/idxstats/{sample}_idxstats.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/idxstats/{sample}_idxstats.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/idxstats/{sample}_idxstats.benchmark"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools view -c -F 4 {input.bam} > {output.counts} 2> >(tee {log.errfile} >&2)
        samtools idxstats {input.bam} > {output.idxstats}  2> >(tee {log.errfile} >&2)
        """


rule assign_virus:
    input:
        idxstats = rules.idxstats.output.idxstats,
        sorted = rules.sort.output.sorted,
        counts = rules.sort.output.counts
    output:
        table = config["inputOutput"]["output_dir"]+"/{sample}/assign_virus/{sample}_count_table.tsv"
    params:
        ref_table = config["resources"]["reference_table"],
        prefix = config["inputOutput"]["output_dir"]+"/{sample}/assign_virus/{sample}_",
        outlier_percentile = config["tools"]["assign_virus"]["outlier_percentile"],
        outlier_percentile_collapsed = config["tools"]["assign_virus"]["outlier_percentile_collapsed"],
        collapse = config["tools"]["assign_virus"]["collapse"],
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/assign_virus/{sample}_assignment.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/assign_virus/{sample}_assignment.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/assign_virus/{sample}_assignment.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/assign_virus.py \
        --ref_table {params.ref_table} \
        --idxstats {input.idxstats} \
        --counts {input.counts} \
        --outprefix {params.prefix} \
        --outlier_percentile {params.outlier_percentile} \
        --outlier_percentile_collapsed {params.outlier_percentile_collapsed} \
        -c {params.collapse}  2> >(tee {log.errfile} >&2)
        """


rule validate_assignment:
    input:
        assignment = rules.assign_virus.output.table,
    output:
        validation = config["inputOutput"]["output_dir"]+"/{sample}/validate_assignment/{sample}_validation.txt",
    params:
        metadata_dir = config["resources"]["metadata_dir"],
        anontable = config["resources"]["anonymization_table"]
    log:
        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/validate_assignment/{sample}_validation.out.log",
        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/validate_assignment/{sample}_validation.err.log",
    benchmark:
        config["inputOutput"]["output_dir"]+"/logs/benchmark/validate_assignment/{sample}_validation.benchmark"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python workflow/scripts/validate_assignment.py \
        --metadata_dir {params.metadata_dir} \
        --anonymization_table {params.anontable} \
        --output {output.validation}  2> >(tee {log.errfile} >&2)
        """
#
#
#rule package_results:
#    input:
#    output:
#    params:
#    log:
#        outfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/package_results/{sample}_package_results.out.log",
#        errfile=config["inputOutput"]["output_dir"]+"/logs/{sample}/package_results/{sample}_package_results.err.log",
#    benchmark:
#        config["inputOutput"]["output_dir"]+"/logs/benchmark/package_results/{sample}_package_results.benchmark"
#    conda:
#        "../envs/python.yaml"
#    shell:
#        "python workflow/scripts/validate_assignment.py --metadata_dir {params.metadata_dir} --anonymization_table {params.anontable} --output {output.validation}"

