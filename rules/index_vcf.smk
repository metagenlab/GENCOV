
rule tbiIndexVCF:
    container: 
        singularity_envs["bcftools"]
    input:
        "{sample}.vcf.gz"
    output:
        temp("{sample}.vcf.gz.tbi")
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
            bcftools index -f -t {input}
        """

rule csiIndexVCF:
    container: 
        singularity_envs["bcftools"]
    input:
        "{sample}.vcf.gz"
    output:
        temp("{sample}.vcf.gz.csi")
    shell:
        r"""
            bcftools index -f {input}
        """

