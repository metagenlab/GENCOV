rule indexPicard:
    input:
        REFERENCE
    output:
        PICARD_INDEX
    conda:
        "../envs/gatk.yaml"
    log:
        os.path.join(PROJFOLDER, "logs", "index", "picard_reference.log")
    shell:
        r"""
            gatk CreateSequenceDictionary -R {input} &> {log}
        """