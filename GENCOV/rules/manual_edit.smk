

rule fasta_edit:
    input:
        fasta_file = os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.iupac_consensus.fasta"),
        vcf = os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.filtered.ALT_corrected.vcf.gz"),
        primers = "data/primers.tsv",
    output:
        # results/consensuses_iupac/freebayes/1015926446_4812.iupac_consensus.edit.fasta
        fasta_edit = os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.iupac_consensus.edit.fasta"),
        modified_pos = os.path.join(IUPAC_CNS_FOLDER, "{snp_calling_tool}", "{sample}.edit_pos.tsv"),
    params:
        frac_filter = CNS_GT_ADJUST,
    script: "../scripts/fasta_edit.py"