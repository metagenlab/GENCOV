singularity: "docker://rkibioinf/samtools:1.11--b05ccf8"

rule getBamStats:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam")
    output:
        stats = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.txt"),
        stats_forR = os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.pipe.txt")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping_stats", "{sample}.bamstats.log")
    conda:
        "../envs/samtools.yaml"
    shell:
        r"""
            samtools flagstat {input.bam} 1> {output.stats} 2> {log};
            cat {output.stats} | sed -e 's/ + /|/' | sed -e 's/ /|/' 1> {output.stats_forR} 2> {log}
        """