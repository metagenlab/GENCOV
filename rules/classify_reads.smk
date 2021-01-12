singularity: "docker://rkibioinf/kraken2:2.1.0--3006bb7"

rule classifyReads:
    input:
        PE1 = os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R1.fastq.gz"),
        PE2 = os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R2.fastq.gz")
    output:
        PE1 = temp(os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R_1.fastq.gz")),
        PE2 = temp(os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R_2.fastq.gz")),
        kraken_out = temp(os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.kraken.out.txt")),
        kraken_report = os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.kraken.report.txt")
    params:
        db = KRAKEN_DB,
        fastq = os.path.join(DATAFOLDER["classified"], "{sample}", "{sample}.R#.fastq.gz")
    log:
        os.path.join(DATAFOLDER["logs"], "classified", "{sample}.classify.log")
    conda:
        "../envs/kraken.yaml"
    threads:
        10
    shell:
        r"""
            ( kraken2 \
                --threads {threads} \
                --db {params.db} \
                --paired \
                --classified-out {params.fastq} \
                --output {output.kraken_out} \
                --report {output.kraken_report} \
                --gzip-compressed \
                {input.PE1} {input.PE2} ) 2> {log}
        """
