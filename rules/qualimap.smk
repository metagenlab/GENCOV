
rule assess_mapping_with_qualimap:
    conda:
        "../envs/qualimap.yaml"   
    input:
        # samples/{sample}/mapping/bwa/{ref}.bam
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam"),
        bam_bai = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam.bai"),
    output:
        report = "report/qualimap/{sample}/qualimapReport.html",
        genome_results = "report/qualimap/{sample}/genome_results.txt",
    log:
        os.path.join(DATAFOLDER["logs"], "{sample}", "qualimap.txt"),
    shell:
        """
        qualimap bamqc -bam {input[bam]} -outdir $(dirname {output[0]}) > {log[0]}
        if [ ! -f {output[report]} ]
        then
            touch {output[report]}
            touch {output[genome_results]}
        fi
        """
