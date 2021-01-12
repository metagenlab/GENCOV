FROM continuumio/miniconda3
ENV VERSION 1.11
ENV TOOL bcftools

# meta data
LABEL base_image="continuumio/miniconda3"
LABEL maintainer="RKI MF1 Bioinformatics <https://www.rki.de/EN/Content/Institute/DepartmentsUnits/MF/MF1/mf1_node.html>"

# install basics
RUN apt update && apt install -y procps wget gzip pigz bc && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels conda-forge && conda config --add channels bioconda && conda config --add channels default

RUN conda install -y -c bioconda $TOOL=$VERSION perl-vcftools-vcf=0.1.16 pbgzip=2016.08.04 && conda clean -a