
def input_trimReads(wildcards):
    
    #if PRIMER:
    #    files = [os.path.join(DATAFOLDER["trimmed"], wildcards.sample, wildcards.sample + ".R1.noprimer.fastq.gz"), 
    #             os.path.join(DATAFOLDER["trimmed"], wildcards.sample, wildcards.sample + ".R2.noprimer.fastq.gz")]
    #else:
    files = list(getFastq(wildcards))
    
    if ADAPTERS:
        files.append(ADAPTERS)
    return files


rule trimReads:
    container: 
        singularity_envs["fastp"]
    input:
        input_trimReads
    output:
        PE1 = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R1.fastq.gz")),
        PE2 = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.R2.fastq.gz")),
        SE1 = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.SE.R1.fastq.gz")),
        SE2 = temp(os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.SE.R2.fastq.gz")),
        json = os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.fastp.json"),
        html = os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.fastp.html")
    log:
        os.path.join(DATAFOLDER["logs"], "trimming", "{sample}.log")
    params:
        adapters = "--adapter_fasta " + ADAPTERS if ADAPTERS else ""
    threads:
        10
    shell:
        r"""
            ( fastp \
                --in1 {input[0]} \
                --out1 {output.PE1} \
                --in2 {input[1]} \
                --out2 {output.PE2} \
                --unpaired1 {output.SE1} \
                --unpaired2 {output.SE2} \
                --json {output.json} \
                --html {output.html} \
                {params.adapters} \
                --qualified_quality_phred 15 \
                --length_required 15 \
                --low_complexity_filter \
                --overrepresentation_analysis \
                --thread {threads} ) &> {log}
        """


