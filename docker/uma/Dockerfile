FROM r-base:4.4.2
ARG buildDate

ADD /UMA_lib/ /UMA_lib/

LABEL about.image="R base image to run UMA panel post-processing for MM alteration detection and analysis"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="r-base:4.4.2"
LABEL image.version="1.0.3"
LABEL buildDate=${buildDate}
LABEL software="Custom R scripts for MM alteration analysis"
LABEL software.version="1.0.4"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \ 
    libcurl4-gnutls-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libssl-dev \
    libssh2-1-dev \
    libxml2-dev \
    procps \
    sudo \
    wget

RUN R -e 'install.packages(c("BiocManager", "data.table", "dplyr", "ggplot2", "gtools", "Hmisc", "mgsub", "NbClust", "optparse", "paletteer", \
		"readr", "readxl", "reshape2", "stringr", "svglite", "systemfonts", "textshaping", "tibble", "tidyr", "vcfR", "XML", "xml2"))'
RUN R -e 'BiocManager::install(c("GenomicRanges", "plyranges"))'

ENTRYPOINT [ "Rscript" ]
