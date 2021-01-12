FROM continuumio/miniconda3
ENV VERSION 1.2.1
ENV TOOL covpipe

# meta data
LABEL base_image="continuumio/miniconda3"
LABEL about.summary="CovPipe is a pipeline to generate consensus sequences from NGS reads based on a reference sequence. The pipeline is tailored to be used for SARS-CoV-2 data, but may be used for other viruses."
LABEL about.license="GLP3"
LABEL about.tags="ncov"
LABEL about.home="https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe"

LABEL maintainer="RKI MF1 Bioinformatics <https://www.rki.de/EN/Content/Institute/DepartmentsUnits/MF/MF1/mf1_node.html>"

# install basics
RUN apt update && apt install -y procps wget gzip pigz bc && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels conda-forge && \
	conda config --add channels bioconda && \
	conda config --add channels default

# download and install covpipe
RUN git clone https://gitlab.com/RKIBioinformaticsPipelines/ncov_minipipe.git
WORKDIR /ncov_minipipe
RUN bash ncov_minipipe.conda.setup/setup_env.sh

# install all environments according to snakemakes default path
RUN for ENV in envs/*.yaml; do MD5=$(md5sum $ENV | awk '{print $1}') && conda env create --prefix .snakemake/conda/$MD5 -f $ENV; done 

# clean conda (keep image small)
RUN conda clean -a

# source environment when creating a container from the image
# docker run --rm covpipe:latest
# or
# docker run --rm covpipe:latest ncov_minipipe --help
WORKDIR /ncov_minipipe
COPY entrypoint.sh .
RUN chmod +x entrypoint.sh
ENTRYPOINT ["/ncov_minipipe/entrypoint.sh"]
#CMD ["ncov_minipipe"]
