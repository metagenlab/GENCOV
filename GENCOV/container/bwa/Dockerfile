FROM continuumio/miniconda3
ENV VERSION 0.7.17
ENV TOOL bwa

# meta data
LABEL base_image="continuumio/miniconda3"
LABEL maintainer="RKI MF1 Bioinformatics <https://www.rki.de/EN/Content/Institute/DepartmentsUnits/MF/MF1/mf1_node.html>"

# install basics
RUN apt update && apt install -y procps wget gzip pigz bc && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels conda-forge && conda config --add channels bioconda && conda config --add channels default

RUN conda install -y -c bioconda $TOOL=$VERSION samtools=1.11 && conda clean -a