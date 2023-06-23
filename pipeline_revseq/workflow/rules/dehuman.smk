## Coutesy of V-pipe: https://github.com/cbg-ethz/V-pipe/tree/master

from functools import partial


rule dh_reuse_alignreject:
    # "rely on aligner's output".
    # this rule re-use the rejected reads in align.smk (e.g. ngshmmalign's /alignments/rejects.sam)
    # (useful when in parallel with the main processing)
    input:
        reject_aln=rules.bwa.output.outfile
    output:
        reject_1=temp("results/{sample}/align_reject/reject_R1.fastq.gz"),
        reject_2=temp("results/{sample}/align_reject/reject_R2.fastq.gz"),
    log:
        outfile="logs/{sample}/align_reject/reject.out.log",
        errfile="logs/{sample}/align_reject/reject.err.log",
    conda:
        config["tools"]["dehuman"]["conda"]
    benchmark:
        "results/{sample}/align_reject/reject.benchmark"
    resources:
        disk_mb=1250,
        #mem_mb=config.bwa_align["mem"],
        #runtime=config.bwa_align["time"],
    threads: config["tools"]["bwa_align"]["threads"]
    shell:
        """
        echo "Keep reject  -----------------------------------------------------"
        echo

        samtools bam2fq -@ {threads} \
                                 -F 2 \
                                 -1 {output.reject_1} \
                                 -2 {output.reject_2} \
                                 {input.reject_aln} \
                                 2> >(tee {log.errfile} >&2)
        echo
        """


rule dh_host_index:
    input:
        reference = config["resources"]["host_ref"]
    output:
        referenceout = config["resources"]["host_ref"] + '.bwt'
    shell:
        "bwa index {input.reference}"


rule dh_hostalign:
    input:
        host_ref=config["resources"]["host_ref"],
        ref_index=rules.dh_host_index.output.referenceout,
        reject_1=rules.dh_reuse_alignreject.output.reject_1,
        reject_2=rules.dh_reuse_alignreject.output.reject_2,
    output:
        host_aln=temp("results/{sample}/dh_hostalign/host_aln.sam"),
    log:
        outfile="logs/{sample}/dh_hostalign/host_aln.out.log",
        errfile="logs/{sample}/dh_hostalign/host_aln.err.log",
    conda:
        config["tools"]["dehuman"]["conda"]
    benchmark:
        "results/{sample}/dh_hostalignhost_aln.benchmark"
    resources:
        disk_mb=1250,
        #mem_mb=config.dehuman["mem"],
        #runtime=config.dehuman["time"],
    threads: config["tools"]["dehuman"]["threads"]
    shell:
        """
        echo "Checking rejects against host's genome  --------------------------"

        bwa mem -t {threads} \
                         -o {output.host_aln}\
                         {input.host_ref} {input.reject_1} {input.reject_2} \
                         > {log.outfile} 2> >(tee {log.errfile} >&2)

        echo
        """


rule dh_filter:
    input:
        host_aln=rules.dh_hostalign.output.host_aln,
        R1=config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_R1.fastq.gz",
        R2=config["inputOutput"]["input_fastqs"] + "/{sample}/{sample}_R2.fastq.gz",
    output:
        filter_count="results/{sample}/dh_filter/dehuman.count",
        filter_list=temp("results/{sample}/dh_filter/dehuman.filter"),
        filtered_1=temp("results/{sample}/dh_filter/filtered_1.fastq.gz"),
        filtered_2=temp("results/{sample}/dh_filter/filtered_2.fastq.gz"),
    params:
        remove_reads_script="../scripts/remove_reads_list.pl",
        keep_host=int(config["tools"]["dehuman"]["keep_host"]),
        sort_tmp=temp("results/{sample}/dh_filter/host_sort.tmp"),
        host_aln_cram="results/{sample}/dh_filter/host_aln.cram",
        F=2,
    log:
        outfile="logs/{sample}/dh_filter/dehuman_filter.out.log",
        errfile="logs/{sample}/dh_filter/dehuman_filter.err.log",
    conda:
        config["tools"]["dehuman"]["conda"]
    benchmark:
        "results/{sample}/dh_filter/dehuman_filter.benchmark"
    resources:
        disk_mb=1250,
        #mem_mb=config.dehuman["mem"],
        #runtime=config.dehuman["time"],
    threads: config["tools"]["dehuman"]["threads"]
    shell:
        """
        # using zcat FILENAME.gz causes issues on Mac, see
        # https://serverfault.com/questions/570024/
        # redirection fixes this:
        unpack_rawreads() {{
            for fq in "${{@}}"; do
                zcat -f < "${{fq}}"
            done
        }}

        echo
        echo "Count aligned reads ---------------------------------------------"
        echo

        count=samtools view -@ {threads} -c -f {params.F} -F 2304 {input.host_aln} | tee {output.filter_count} 2> >(tee {log.errfile} >&2) )

        if (( count > 0 )); then
            echo
            echo "-----------------------------------------------------------------"
            echo "Needs special care: ${{count}} potential human reads found"
            echo "-----------------------------------------------------------------"
            echo
            echo "Removing identified host reads from raw reads -------------------"
            echo

            # get list
            samtools view -@ {threads} \
                                   -f {params.F} \
                                   {input.host_aln} \
                                   | cut -f 1 > {output.filter_list} 2> >(tee -a {log.errfile} >&2)

            unpack_rawreads {input.R1:q} \
                   | {params.remove_reads_script} {output.filter_list} \
                   | gzip \
                   > {output.filtered_1} 2> >(tee -a {log.errfile} >&2) &

            unpack_rawreads {input.R2:q} \
                   | {params.remove_reads_script} {output.filter_list} \
                   | gzip \
                   > {output.filtered_2} 2> >(tee -a {log.errfile} >&2) &

            wait

            if (( {params.keep_host} )); then
                # keep the rejects for further analysis

                echo
                echo "Keeping host-aligned virus' rejects ------------------------------"
                echo


                # (we compress reference-less, because the reference size is larger
                # than the contaminant reads)
                rm -f '{params.sort_tmp}'.[0-9]*.bam
                FMT=cram,no_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000
                samtools view -@ {threads} \
                                      -h -f 2 \
                                      {input.host_aln} \
                    | samtools sort -@ {threads} \
                                    -T {params.sort_tmp} \
                                    --output-fmt ${{FMT}} \
                                    -o {params.host_aln_cram} \
                                    2> >(tee -a {log.errfile} >&2)

                echo
                echo "Compressing host-depleted raw reads ------------------------------"
                echo

                samtools index -@ {threads} {params.host_aln_cram} 2> >(tee -a {log.errfile} >&2)
            fi
        else
            echo
            echo "No potential human reads found -----------------------------------"
            echo "Copy raw reads file"
            echo
            unpack_rawreads {input.R1} | gzip > {output.filtered_1} 2> >(tee -a {log.errfile} >&2) &
            unpack_rawreads {input.R2} | gzip > {output.filtered_2} 2> >(tee -a {log.errfile} >&2) &
            wait
            touch {output.filter_list}
        fi
        echo
        """


rule dehuman:
    input:
        global_ref=config["resources"]["reference"],
        ref_index=config["resources"]["reference"] + '.bwt',
        filtered_1=rules.dh_filter.output.filtered_1,  # =temp_prefix("{dataset}/raw_uploads/filtered_1.fastq.gz"),
        filtered_2=rules.dh_filter.output.filtered_2,  # =temp_prefix("{dataset}/raw_uploads/filtered_2.fastq.gz"),
    output:
        cram_sam=temp("results/{sample}/dehuman/dehuman.sam"),
        final_cram="results/{sample}/dehuman/dehuman.cram",
        checksum="results/{sample}/dehuman/dehuman.cram.%s" % config["tools"]["dehuman"]["checksum"],
    params:
        checksum_type=config["tools"]["dehuman"]["checksum"],
        sort_tmp=temp("results/{sample}/dehuman/dehuman.tmp"),
    log:
        outfile="logs/{sample}/dehuman/dehuman.out.log",
        errfile="logs/{sample}/dehuman/dehuman.err.log",
    conda:
        config["tools"]["dehuman"]["conda"]
    benchmark:
        "results/{sample}/dehuman/dehuman.benchmark"
    resources:
        disk_mb=1250,
        #mem_mb=config.dehuman["mem"],
        #runtime=config.dehuman["time"],
    threads: config["tools"]["dehuman"]["threads"]
    shell:
        """
        echo "Compress filtered sequences --------------------------------------"
        echo

        bwa mem -t {threads} \
                         -C \
                         -o {output.cram_sam} \
                         {input.global_ref} {input.filtered_1} {input.filtered_2} \
                         > {log.outfile} 2> >(tee {log.errfile} >&2)

        # HACK handle incompatibilities between:
        #  - Illumina's 'bcl2fastq', which write arbitrary strings
        #  - 'bwa mem' which keep comments in the SAM file verbatim as in the FASTQ file
        #  - 'samtools' which expects comment to be properly marked as 'BC:Z:'
        #    as per SAM format specs
        REGEXP=\'s{{(?<=\\t)([[:digit:]]:[[:upper:]]:[[:digit:]]:([ATCGN]+(\+[ATCGN]+)?|[[:digit:]]+))$}}{{BC:Z:\\1}}\'
        FMT=cram,embed_ref,use_bzip2,use_lzma,level=9,seqs_per_slice=1000000

        rm -f '{params.sort_tmp}'.[0-9]*.bam
        perl -p -e ${{REGEXP}} {output.cram_sam} \
              | samtools sort -@ {threads} \
                                       -T {params.sort_tmp} \
                                       -M \
                                       --reference {input.global_ref} \
                                       --output-fmt ${{FMT}} \
                                       -o {output.final_cram} \
                                       2> >(tee -a {log.errfile} >&2)

        {params.checksum_type}sum {output.final_cram} > {output.checksum} 2> >(tee -a {log.errfile} >&2)

        echo
        echo DONE -------------------------------------------------------------
        echo
        """
