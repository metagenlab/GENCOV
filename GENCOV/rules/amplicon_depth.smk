
rule amplicon_depth:
    container: 
        singularity_envs["bedtools"]
    input:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_{sample}.bam")
    output:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "amplicon_depth.{sample}.tsv")
    params:
        bed_file = config["amplicons_bed"]
    shell:
        """
        # get depth at each position
        bedtools coverage -b {input[0]} -a {params[0]} -d > {output[0]}
        """

rule stats_amplicon_depth:
    input:
        expand(os.path.join(DATAFOLDER["mapping"], "{sample}", "amplicon_depth.{sample}.tsv"), sample=SAMPLES)
    output: 
        "results/amplicon_depth.tsv",
        "results/samples_depth.tsv"
    script:
        "../scripts/stats_amplicon_depth.py"