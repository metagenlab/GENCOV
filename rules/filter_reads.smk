singularity: "docker://rkibioinf/general:3.6.0--28150df"

rule getVirusReads:
    input:
        PE1 = os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R_1.fastq.gz"),
        PE2 = os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R_2.fastq.gz")
    output:
        PE1 = temp(os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R1.fastq.gz")),
        PE2 = temp(os.path.join(DATAFOLDER["filtered"], "{sample}", "{sample}.R2.fastq.gz"))
    params:
        taxid = KRAKEN_TAX_ID
    log:
        os.path.join(DATAFOLDER["logs"], "filtered", "{sample}.extract.log")
    shell:
        r"""
            set +o pipefail
            zgrep -A3 'kraken:taxid|{params.taxid}' {input.PE1} | sed -e 's/^--$//' | sed '/^\s*$/d' | gzip 1> {output.PE1} 2>> {log}
            zgrep -A3 'kraken:taxid|{params.taxid}' {input.PE2} | sed -e 's/^--$//' | sed '/^\s*$/d' | gzip 1> {output.PE2} 2>> {log}
        """
