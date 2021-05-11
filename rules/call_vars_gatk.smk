

rule callGenomicVariants_gatk:
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_small_align_{sample}.sort.bam"),
        ref = REFERENCE,
        ref_dict = re.sub(".fasta", "",REFERENCE) + '.dict',
        ref_index = REFERENCE + ".fai",
        bam_index = os.path.join(DATAFOLDER["mapping"], "{sample}", "filtered_primers_small_align_{sample}.sort.bam.bai"),
    output:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "gatk", "{sample}.vcf")
    params:
        cov = VAR_CALL_COV,
        vcount = VAR_CALL_COUNT,
        frac = VAR_CALL_FRAC,
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "gatk", "{sample}.log")
    conda:
        "../envs/gatk.yaml"
    threads: 1
    shell:
        r"""
            #freebayes -f {input.ref} --min-alternate-count {params.vcount} --min-alternate-fraction {params.frac} --min-coverage {params.cov} --pooled-continuous --haplotype-length -1 {input.bam} 1> {output.vcf} 2> {log}
            gatk HaplotypeCaller --max-reads-per-alignment-start 0 --sample-ploidy 1 --output {output.vcf} --input {input.bam} --reference {input.ref} 2> {log}
        """

