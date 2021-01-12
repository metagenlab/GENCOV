singularity: "docker://rkibioinf/bcftools:1.11--19c96f3"

rule tbiIndexVCF:
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
    input:
        "{sample}.vcf.gz"
    output:
        temp("{sample}.vcf.gz.csi")
    conda:
        "../envs/bcftools.yaml"
    shell:
        r"""
            bcftools index -f {input}
        """

