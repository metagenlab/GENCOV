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
RUN git clone https://github.com/metagenlab/GENCOV.git && echo OK
WORKDIR /GENCOV
RUN bash ncov_minipipe.conda.setup/setup_env.sh && for ENV in envs/*.yaml; do MD5=$(md5sum $ENV | awk '{print $1}') && conda env create --prefix .snakemake/conda/$MD5 -f $ENV; done  &&  conda clean -a

ENTRYPOINT ["/bin/bash"]
ENV PATH /GENCOV/covpipe_environment/bin:$PATH