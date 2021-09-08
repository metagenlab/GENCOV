
rule indexSamtools:
    container: 
        singularity_envs["samtools"]
    input:
        "{fname}"
    output:
        "{fname}.fai"
    shell:
        r"""
            samtools faidx {input}
        """
