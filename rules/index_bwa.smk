
rule indexBwa:
    container: 
        singularity_envs["bwa"]
    input:
        "{fname}"
    output:
        "{fname}.bwt"
    shell:
        r"""
            bwa index {input} &> /dev/null
        """