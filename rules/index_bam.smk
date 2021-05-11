singularity: "docker://rkibioinf/samtools:1.11--b05ccf8"

rule indexMapping:
    input:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam")
    output:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam.bai")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping", "index_{sample}.log")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            samtools index {input} &> {log}
        """