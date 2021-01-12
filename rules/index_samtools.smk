singularity: "docker://rkibioinf/samtools:1.11--b05ccf8"

rule indexSamtools:
    input:
        "{fname}"
    output:
        "{fname}.fai"
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            samtools faidx {input}
        """
