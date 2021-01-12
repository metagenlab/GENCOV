singularity: "docker://rkibioinf/snpeff:4.5covid19--c40241f"

rule annotateVariants:
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.vcf")
    output:
        report = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.annotation.html"),
        ann_vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{sample}.annotation.covered.af.vcf")
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.annotate.vcf.log")
    params:
        reference = REFERENCE
    conda:
        "../envs/snpeff.yaml"
    threads:
        1
    shell:
        r"""
            # get genome name
            genome_name=$(head -n1 {params.reference} | cut -f1 -d' ' | sed 's/>//')

            # run annotation
            # -t {threads} \
            snpEff ann \
                -noLog \
                -stats {output.report} \
                $genome_name \
                {input} 1> {output.ann_vcf} 2> {log}
        """
