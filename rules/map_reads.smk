singularity: "docker://rkibioinf/bwa:0.7.17--a9f152d"

def input_map2reference(wildcards):
    if KRAKEN_DB:
        PE1 = os.path.join(DATAFOLDER["filtered"], wildcards.sample, wildcards.sample + ".R1.fastq.gz")
        PE2 = os.path.join(DATAFOLDER["filtered"], wildcards.sample, wildcards.sample + ".R2.fastq.gz")
    else:
        PE1 = os.path.join(DATAFOLDER["trimmed"], wildcards.sample, wildcards.sample + ".R1.fastq.gz")
        PE2 = os.path.join(DATAFOLDER["trimmed"], wildcards.sample, wildcards.sample + ".R2.fastq.gz")
    return {
            'PE1': PE1, 
            'PE2': PE2 
           }

rule map2reference:
    # add read group tags: https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups
    input:
        unpack(input_map2reference),
        ref = REFERENCE,
        ref_index = BWA_INDEX
    output:
        os.path.join(DATAFOLDER["mapping"], "{sample}", "unfiltered_{sample}.bam")
    log:
        os.path.join(DATAFOLDER["logs"], "mapping", "{sample}.log")
    conda:
        "../envs/bwa.yaml"
    threads:
        10
    shell:
        r"""
            (   time \
                bwa mem -t {threads} \
                    -R '@RG\tID:{wildcards.sample}\tPU:{wildcards.sample}\tSM:{wildcards.sample}\tPL:ILLUMINA\tLB:000' \
                    {input.ref} \
                    {input.PE1} {input.PE2} | \
                    samtools view -Sb -@ {threads} -o {output} \
            ) &> {log}
        """


