singularity: "docker://rkibioinf/samtools:1.11--b05ccf8"

rule sortBam:
    input:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_small_align_{sample}.bam")
    output:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_small_align_{sample}.sort.bam")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping", "{sample}.sort.log")
    conda:
        "../envs/samtools.yaml"
    threads:
        10
    shell:
        r"""
            samtools sort -@ {threads} -o {output} {input} &> {log}
        """