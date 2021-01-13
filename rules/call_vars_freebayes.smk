
singularity: "docker://rkibioinf/freebayes:1.3.2--1793b52"

rule callGenomicVariants_freebayes:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam"),
        ref = REFERENCE,
        ref_index = REFERENCE + ".fai"
    output:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "freebayes", "{sample}.vcf")
    params:
        cov = VAR_CALL_COV,
        vcount = VAR_CALL_COUNT,
        frac = VAR_CALL_FRAC,
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "freebayes", "{sample}.log")
    conda:
        "../envs/freebayes.yaml"
    threads: 1
    shell:
        r"""
        freebayes -f {input.ref} --min-alternate-count {params.vcount} --min-alternate-fraction {params.frac} --min-coverage {params.cov} --pooled-continuous --haplotype-length -1 {input.bam} 1> {output.vcf} 2> {log}
        """

