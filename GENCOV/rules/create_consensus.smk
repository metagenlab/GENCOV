# low coverage hard filtering
rule createMaskConsensus:
    container: 
        singularity_envs["bedtools"]
    input:
        bam = os.path.join(DATAFOLDER["mapping"], "{sample}", "{sample}.sort.bam")
    output:
        os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.lowcov.raw.bed")
    params:
        cov = CNS_MIN_COV
    log:
        os.path.join(DATAFOLDER["logs"], "masking", "{sample}.consensus_mask.log")
    shell:
        r"""
            bedtools genomecov -bga -ibam {input.bam} | awk '$4 < {params.cov}' | bedtools merge  > {output[0]}
        """

# retrieve deletion positions
rule extract_del_positions:
    container: 
        singularity_envs["bcftools"]
    input:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.vcf")
    output:
        os.path.join(DATAFOLDER["masking"], "{sample}", "{snp_calling_tool}","{sample}.del.bed")
    params:
        cov = CNS_MIN_COV
    log:
        os.path.join(DATAFOLDER["logs"], "masking", "{sample}.{snp_calling_tool}.consensus_mask.log")
    shell:
        r"""
            bcftools view -v indels {input[0]} | bcftools view -H | awk '{{print $1,$2,$2+length($4)}}' | sed 's/ /\t/g' > {output[0]}
        """

# do not filter deletion positions
rule filter_del_positions:
    container: 
        singularity_envs["bedtools"]
    input:
        os.path.join(DATAFOLDER["masking"], "{sample}", "{sample}.lowcov.raw.bed"),
        os.path.join(DATAFOLDER["masking"], "{sample}", "{snp_calling_tool}", "{sample}.del.bed")
    output:
        os.path.join(DATAFOLDER["masking"], "{sample}", "{snp_calling_tool}", "{sample}.lowcov.bed")
    params:
        cov = CNS_MIN_COV
    log:
        os.path.join(DATAFOLDER["logs"], "masking", "{sample}.{snp_calling_tool}.consensus_mask.log")
    shell:
        r"""
            bedtools subtract -a {input[0]} -b {input[1]} > {output[0]}
        """


## var hard filtering
rule filterVarsConsensus_gatk:
    container: 
        singularity_envs["bcftools"]
    input:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "gatk", "{sample}.vcf")
    output:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "gatk", "{sample}.filtered.vcf.gz")
    params:
        mqm = VAR_FILTER_MQM,
        sap = VAR_FILTER_SAP,
        qual = VAR_FILTER_QUAL
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.gatk.filtered.vcf.log")
    shell:
        r"""
			temp_out="{output}"
			temp_out="${{temp_out%%.vcf.gz}}.vcf"
			{{
				set -x 
				bcftools filter -e \
					"INFO/MQ < {params.mqm} | INFO/FS > {params.sap} | QUAL < {params.qual}" \
					 -o "{output}" -O z "{input}"
			}} |& tee {log} 
        """


## var hard filtering
rule filterVarsConsensus_freebayes:
    container: 
        singularity_envs["bcftools"]
    input:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "freebayes", "{sample}.vcf")
    output:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "freebayes", "{sample}.filtered.vcf.gz")
    params:
        mqm = VAR_FILTER_MQM,
        sap = VAR_FILTER_SAP,
        qual = VAR_FILTER_QUAL
    log:
        os.path.join(DATAFOLDER["logs"], "variant_calling", "{sample}.freebayes.filtered.vcf.log")
    shell:
        r"""
			temp_out="{output}"
			temp_out="${{temp_out%%.vcf.gz}}.vcf"
			{{
				set -x 
				bcftools filter -e \
					"INFO/MQM < {params.mqm} | INFO/SAP > {params.sap} | QUAL < {params.qual}" \
					 -o "{output}" -O z "{input}"
			}} |& tee {log} 
        """

def get_ad_tag(wildcards):
    if wildcards.snp_calling_tool == 'gatk':
        return 'AD'
    if wildcards.snp_calling_tool == 'freebayes':
        return 'AD'
    else:
        raise IOError("Unknown snp caller")



## genotype adjustment 
rule adjustGtConsensus:
    container: 
        singularity_envs["pyvcf"]
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.vcf.gz")
    output:
        os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.ALT_corrected.vcf"),
    params:
        ao_tag = get_ad_tag,
        dp_tag = "DP",
        min_freq=0.1,
    script:
        "../scripts/update_ALT_freq.py"


## filter fariants based on ALT frequency 
rule filter_ALTF_freq:
    input:
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.ALT_corrected.vcf.gz"),
    output:
        filtered_vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.ALT_corrected.freq_filter.vcf"),
        modified_pos = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.edit_pos.tsv"),
    params:
        freq = CNS_GT_ADJUST,
        frac_filter = CNS_GT_ADJUST,
        keep_positions = config["mutation_filter"]
    script: "../scripts/filter_vcf.py"


## filter fariants based on ALT frequency 
rule bgzip_vcf:
    container: 
        singularity_envs["bgzip"]
    input:
        "{path}.vcf",
    output:
        "{path}.vcf.gz",
    shell:
        """
        bgzip -c {input[0]} > {output[0]}
        """

## create ambig consensus
## use filtered vcf with alt freq > CNS_GT_ADJUST
def input_createAmbiguousConsensus(wildcards):
    files = {}
    files['fasta'] = REFERENCE
    files['mask'] = os.path.join(DATAFOLDER["masking"], wildcards.sample, "{snp_calling_tool}", wildcards.sample + ".lowcov.bed")
    files['vcf'] = os.path.join(DATAFOLDER["variant_calling"], wildcards.sample, "{snp_calling_tool}", wildcards.sample + ".filtered.ALT_corrected.freq_filter.vcf.gz")
    files['vcf_index'] = os.path.join(DATAFOLDER["variant_calling"], wildcards.sample, "{snp_calling_tool}", wildcards.sample + ".filtered.ALT_corrected.freq_filter.vcf.gz.tbi")
    return files

rule createAmbiguousConsensus:
    container: 
        singularity_envs["bcftools"]
    input:
        unpack(input_createAmbiguousConsensus)
    output:
        temp(os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.iupac_consensus.tmp"))
    log:
        os.path.join(DATAFOLDER["logs"], "consensus", "{sample}.{snp_calling_tool}.iupac_consensus.tmp.log")
    shell:
        r"""
            
            ( bcftools consensus \
                -m {input.mask} \
                -o {output} \
                -f {input.fasta} \
                --sample {wildcards.sample} \
                {input.vcf} ) &> {log}
        """

## adapt consensus header
rule createHeaderConsensus:
    input:
        fasta = os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.iupac_consensus.tmp"),
    output:
        os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.iupac_consensus.fasta")
    log:
        os.path.join(DATAFOLDER["logs"], "consensus", "{sample}.{snp_calling_tool}.iupac_consensus.fasta.log")
    shell:
        r"""
            sed "1 s/.*/>{wildcards.sample}/" {input.fasta} 1> {output} 2> {log}
        """

## create masked consensus
rule createMaskedConsensus:
    input:
        fasta = os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.iupac_consensus.tmp"),
        version = os.path.join(PROJFOLDER, "pipeline.version")
    output:
        os.path.join(MASKED_CNS_FOLDER, "{snp_calling_tool}", "{sample}.masked_consensus.fasta")
    log:
        os.path.join(DATAFOLDER["logs"], "consensus", "{sample}.{snp_calling_tool}.masked_consensus.log")
    shell:
        r"""
            echo ">{wildcards.sample}_masked_consensus\n" 1> {output} 2> {log}
            tail -n +2 {input.fasta} | tr "RYSWKMBDHVN" "N" 1>> {output} 2>> {log}
        """
