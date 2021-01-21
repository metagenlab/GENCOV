


rule countmixed:
    input:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.mixed.vcf"),
    output:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.mixed_count.txt"),
    conda:
        "../envs/bcftools.yaml"
    singularity: 
        "docker://rkibioinf/bcftools:1.11--19c96f3"
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "{snp_calling_tool}.{sample}.gatk.filtered.mixed.log")
    shell:
        r"""
			{{
                n_mixed=`bcftools query -f '%POS\t%REF\t%ALT\t[%DP\t%RO\t%AO\t%AF]\n' {input[0]} | wc -l`;
                echo -e {wildcards.sample}\\t$n_mixed > {output}
			}} |& tee {log} 
        """


rule combine_mixed_counts:
    input:
        expand(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{{snp_calling_tool}}", "{sample}.mixed_count.txt"), sample=SAMPLES),
    output:
        "results/{snp_calling_tool}_mixed_count.tsv",
    shell:
        r"""
			cat {input} >> {output}
        """