# cd ${params.sd} command problematic in container

rule getVersion:
    output:
        os.path.join(PROJFOLDER, "pipeline.version")
    params:
        sd = srcdir(".")
    log:
        os.path.join(DATAFOLDER["logs"], "version", "git.repository.version.log")
    shell:
        r"""
            echo git version
            wd=$PWD
            which git
            cd {params.sd};
            git describe 1> $wd/{output} 2> $wd/{log} || echo 'unknown_version' 1> $wd/{output} 2> $wd/{log}
        """