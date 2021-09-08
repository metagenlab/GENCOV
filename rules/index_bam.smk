
rule indexMapping:
    container: 
        singularity_envs["samtools"]
    input:
        "{any}.bam"
    output:
        "{any}.bam.bai"
    log:
        os.path.join(DATAFOLDER["logs"], "mapping", "{any}.index.log")
    shell:
        r"""
            samtools index {input} &> {log}
        """