rule prepareReference:
    input:
        ref = config['reference']
    output:
        ref = REFERENCE
    shell:
        r"""
            mkdir -p $(dirname "{output.ref}")
            cat {input.ref} > {output.ref}  
            #dos2unix {output.ref}
        """

