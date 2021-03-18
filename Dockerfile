FROM continuumio/miniconda3:4.7.12

# meta data
LABEL base_image="continuumio/miniconda3:4.7.12"
LABEL about.tags="GENCOV"

LABEL maintainer="Trestan Pillonel"

# install basics
RUN apt update && apt install -y procps wget gzip pigz bc && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels conda-forge && \
	conda config --add channels bioconda && \
	conda config --add channels default

# download and install covpipe
RUN git clone https://github.com/metagenlab/GENCOV.git && echo OK6
WORKDIR /GENCOV
RUN bash ncov_minipipe.conda.setup/setup_env.sh
RUN /GENCOV/covpipe_environment/bin/snakemake -s /GENCOV/ncov_minipipe.snake --cores 16 --use-conda --conda-prefix /GENCOV/.snakemake/conda --configfile /GENCOV/ncov_minipipe.config --conda-create-envs-only  &&  conda clean -a

# --conda-create-envs-only --conda-prefix ${conda_folder} epidemiology --config species="Mycobacterium_tuberculosis"


ENTRYPOINT ["/bin/bash"]
ENV PATH /GENCOV/covpipe_environment/bin:$PATH