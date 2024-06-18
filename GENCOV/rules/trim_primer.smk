
if config["trimming"] == 'fgbio':

    rule remove_primers:
        container: 
            singularity_envs["fgbio"]
        # add read group tags: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "unfiltered_{sample}.bam"),
            REFERENCE,
            config["primer"]
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_{sample}.bam")
        log:
            os.path.join(DATAFOLDER["logs"], "mapping", "fgbio_{sample}.log")
        threads:
            10
        shell:
            r"""
                (   time \
                    fgbio \
                        --sam-validation-stringency=LENIENT \
                        TrimPrimers \
                        -i {input[0]} \
                        -o {output[0]} \
                        -p {input[2]} \
                        -r {input[1]} \
                        -H true
                ) &> {log}
            """


if config["trimming"] == 'ivar':

    rule remove_primers:
        container: 
            singularity_envs["ivar"]
        # add read group tags: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "unfiltered_{sample}.bam"),
            REFERENCE,
            config["primer"],
            os.path.join(DATAFOLDER["mapping"], "{sample}", "unfiltered_{sample}.bam.bai"),
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_{sample}.ivar.bam")
        log:
            os.path.join(DATAFOLDER["logs"], "mapping", "fgbio_{sample}.log")
        threads:
            10
        shell:
            r"""
                (   time \
                    ivar trim -i {input[0]} -b {input[2]} -m 20 -e -p {output[0]}
                ) &> {log}
            """

    rule sort_ivar:
        container: 
            singularity_envs["samtools"]
        input:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_{sample}.ivar.bam")
        output:
            os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_{sample}.bam")
        log:
            os.path.join(DATAFOLDER["logs"], "mapping", "{sample}.ivar.sort.log")
        threads:
            10
        shell:
            r"""
                samtools sort -@ {threads} -o {output} {input} &> {log}
            """


rule filter_bam_small_alignments:
    container: 
        singularity_envs["fgbio"]
    """
    filter alignments smaller than 20bp from bam file
    """
    input:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_{sample}.bam")
    output:
        temp(os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_small_align_{sample}.bam"))
    shell:  
        """
        fgbio FilterBam -i {input[0]} -m 20 -o {output[0]}
        """
