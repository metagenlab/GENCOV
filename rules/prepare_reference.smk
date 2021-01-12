singularity: "docker://rkibioinf/dos2unix:7.4.1--34af548"

rule prepareReference:
    input:
        ref = config['reference']
    output:
        ref = REFERENCE
    conda:
        "../envs/dos2unix.yaml"
    shell:
        r"""
            mkdir -p $(dirname "{output.ref}")
            cat {input.ref} > {output.ref}  
            dos2unix {output.ref}
        """

