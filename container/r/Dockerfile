FROM continuumio/miniconda3
ENV VERSION 0.1
ENV TOOL rscript

# meta data
LABEL base_image="continuumio/miniconda3"
LABEL maintainer="RKI MF1 Bioinformatics <https://www.rki.de/EN/Content/Institute/DepartmentsUnits/MF/MF1/mf1_node.html>"

# install basics
RUN apt update && apt install -y procps wget gzip pigz bc && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# configure conda channels
RUN conda config --add channels conda-forge && conda config --add channels bioconda && conda config --add channels r && conda config --add channels default

RUN conda install -y pandoc=2.2.3.2 pango=1.42.4 freetype=2.9.1 tk=8.6.8 harfbuzz=1.8.8 r=3.6.0 \
pixman=0.38 zstd=1.3.7 gmp=6.1.2 r-data.table=1.12.2 r-formattable=0.2.0.1 r-ggplot2=3.3.2 \
r-kableextra=1.3.1 r-knitr=1.30 r-reshape2=1.4.3 r-rjson=0.2.20 r-rmarkdown=2.5 mscorefonts=0.0.1 \
&& conda clean -a
