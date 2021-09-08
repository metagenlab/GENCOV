
rule getInsertSize:
    container: 
        singularity_envs["samtools"]
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam")
    output:
        sizes = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.fragsize.tsv")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping_stats", "{sample}.fragsize.log")
    shell:
        r"""
            echo 0 1> {output.sizes}
            samtools view -F 4 {input.bam} | cut -f9 1>> {output.sizes} 2> {log}
        """