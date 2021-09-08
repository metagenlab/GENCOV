
rule indexSamtools:
    container: 
        singularity_envs["samtools"]
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
