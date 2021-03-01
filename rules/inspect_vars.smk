singularity: "docker://rkibioinf/snpeff:4.5covid19--c40241f"

rule annotateVariants:
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.ALT_corrected.vcf.gz"),
    output:
        report = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.annotation.html"),
        csv = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.stats.csv"),
        ann_vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.annotation.covered.af.vcf")
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.{snp_calling_tool}.annotate.vcf.log")
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
                -hgvs1LetterAa \
                -stats {output.report} \
                -csvStats {output.csv} \
                $genome_name \
                {input} 1> {output.ann_vcf} 2> {log}
        """

rule VariantTable:
    input:
        ann_vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.annotation.covered.af.vcf")
    output:
        snp_table = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.annotation.tsv")
    conda:
        "../envs/snpsift.yaml"
    shell:
        """
        SnpSift extractFields {input[0]} CHROM POS REF ALT 'ANN[0].EFFECT' 'ANN[0].GENE' 'ANN[0].HGVS_P' 'ANN[0].AA_POS' DP RO AO > {output[0]}
        """


rule merge_variant_tables:
    input:
        all_tables = expand(os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{{snp_calling_tool}}", "{sample}.annotation.tsv"), sample=SAMPLES)
    output:
        "results/{snp_calling_tool}_all_variants.tsv",
        "results/{snp_calling_tool}_all_variants_summary.tsv"
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/combine_tables.py"
        