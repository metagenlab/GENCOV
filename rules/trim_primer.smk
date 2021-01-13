singularity: "docker://rkibioinf/ptrimmer:1.3.3--a577b16"

def input_trimPrimer(wildcards):
    files = list(getFastq(wildcards))
    if PRIMER:
        files.append(PRIMER)
    return files

rule trimPrimer:
    input:
        input_trimPrimer
    output:
        PE1 = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R1.noprimer.fastq")),
        PE2 = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R2.noprimer.fastq")),
        PE1gz = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R1.noprimer.fastq.gz")),
        PE2gz = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R2.noprimer.fastq.gz"))
    log:
        os.path.join(DATAFOLDER["logs"], "trimming", "{sample}.primerremoval.log")
    params:
        m = PRIMER_MISMATCH
    conda:
        "../envs/ptrimmer.yaml"
    threads: 1
    shell:
        r"""
            ptrimmer \
                -t pair \
                -a {input[2]} \
                -f {input[0]} \
                -d {output.PE1} \
                -r {input[1]} \
                -e {output.PE2} \
                -m {params.m} \
                --keep &> {log}
            gzip -k {output.PE1} {output.PE2}
        """