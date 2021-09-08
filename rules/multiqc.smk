
rule copy_result_files_for_multiqc:
    input:
        fastp_json = os.path.join(DATAFOLDER["trimmed"], "{sample}", "{sample}.fastp.json"),
        qualimap_txt =  "report/qualimap/{sample}/genome_results.txt",
        qualimap_report =  "report/qualimap/{sample}/qualimapReport.html",
        snpeff_report=os.path.join(DATAFOLDER["variant_calling"], "{sample}", "{snp_calling_tool}", "{sample}.stats.csv"),
    output:
        "report/multiqc/{snp_calling_tool}/log/{sample}.log.txt"
    shell:
        """
        mkdir -p $(dirname {output[0]})/{wildcards.sample}
        cp {input[0]} $(dirname {output[0]})/{wildcards.sample}
        cat {input[1]} | sed 's/.sort.bam/.bam/g' > $(dirname {output[0]})/{wildcards.sample}/genome_results.txt
        cat {input[2]} | sed 's/.sort.bam/.bam/g' > $(dirname {output[0]})/{wildcards.sample}/qualimapReport.html
        mkdir $(dirname {output[0]})/{wildcards.sample}/raw_data_qualimapReport
        for onefile in `ls $(dirname {input[1]})/raw_data_qualimapReport/*txt`; do 
            cat $onefile > $(dirname {output[0]})/{wildcards.sample}/raw_data_qualimapReport/$(basename $onefile);
        done
        cp {input[3]} $(dirname {output[0]})/{wildcards.sample}
        echo "okay" > {output}
        """


rule multiqc:
    container: 
        singularity_envs["multiqc"]
    input:
        expand("report/multiqc/{{snp_calling_tool}}/log/{sample}.log.txt", sample=SAMPLES)
    output:
        "report/multiqc/{snp_calling_tool}/multiqc_report.html",
        "report/multiqc/{snp_calling_tool}/multiqc_data/multiqc_data.json"
    log:
        os.path.join(DATAFOLDER["logs"], "reporting", "{snp_calling_tool}.multiqc.log")
    threads:
        1
    shell:
        r"""
        multiqc -f -o $(dirname {output[0]}) $(dirname {input[0]} | tr "\n" " ") &> {log[0]}
        """
