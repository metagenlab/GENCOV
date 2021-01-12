singularity: "docker://rkibioinf/bwa:0.7.17--a9f152d"

rule indexBwa:
    input:
        "{fname}"
    output:
        "{fname}.bwt"
    conda:
        "../envs/bwa.yaml"
    shell:
        r"""
            bwa index {input} &> /dev/null
        """