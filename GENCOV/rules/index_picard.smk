rule indexPicard:
    container: 
        singularity_envs["gatk"]
    input:
        REFERENCE
    output:
        PICARD_INDEX
    log:
        os.path.join(PROJFOLDER, "logs", "index", "picard_reference.log")
    shell:
        r"""
            gatk CreateSequenceDictionary -R {input} &> {log}
        """