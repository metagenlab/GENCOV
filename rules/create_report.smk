singularity: "docker://rkibioinf/rscript:0.1--71f6226"

rule createReport:
    input:
        coverage = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.coverage.tsv"), sample = SAMPLES.keys()),
        frag_size  = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.fragsize.tsv"), sample = SAMPLES.keys()),
        mapping_statistics = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.txt"), sample = SAMPLES.keys()),
        mapping_statistics_forR = expand(os.path.join(DATAFOLDER["mapping_stats"], "{sample}", "{sample}.bamstats.pipe.txt"), sample = SAMPLES.keys()),
        version = os.path.join(PROJFOLDER, "pipeline.version")
    output:
        report = os.path.join(PROJFOLDER, "qc_report.html"),
        csv = os.path.join(DATAFOLDER["reporting"], "coverage_samples.csv")
    params:
        p_folder = PROJFOLDER,
        l_folder = DATAFOLDER["reporting"],
        run_id = REPORT_RUNID,
        tax_id = KRAKEN_TAX_ID,
        template = srcdir("../ncov_minipipe.Rmd")
    log:
        os.path.join(DATAFOLDER["logs"], "reporting", "reporting.log")
    conda:
        "../envs/r.yaml"
    threads:
        1
    shell:
        # maybe need to replace shell by r call
        r"""
            # create report
            echo "####### compiling report" >> {log}
            VERSION=$(cat {input.version})
            Rscript -e "rmarkdown::render('{params.template}',
                                            params=list(proj_folder='{params.p_folder}', list_folder='{params.l_folder}', run_name='{params.run_id}', tax_id='{params.tax_id}', version='$VERSION'),
                                            output_file=file.path('{output.report}'))" &> {log}
        """
