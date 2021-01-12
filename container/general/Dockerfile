FROM continuumio/miniconda3
ENV VERSION 3.6.0
ENV TOOL general

# meta data
LABEL base_image="continuumio/miniconda3"
LABEL maintainer="RKI MF1 Bioinformatics <https://www.rki.de/EN/Content/Institute/DepartmentsUnits/MF/MF1/mf1_node.html>"

# install basics
RUN apt update && apt install -y procps wget gzip pigz bc git && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels bioconda && conda config --add channels anaconda && conda config --add channels conda-forge &&  conda config --add channels default

RUN conda install -y python=$VERSION snakemake=5.3.0 strictyaml=1.1.0 biopython=1.78 pandas=1.1.4 && conda clean -a