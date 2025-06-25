FROM r-base:4.4.2
ARG buildDate

ADD /ichorCNA-lib/ /ichorCNA/

LABEL about.image="R base image to run ichorCNA for CN alterations estimation on ULP bam files"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="r-base:4.4.2"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="ichorCNA for ULP-WGS cfDNA CN alteration detection and tumor fraction estimation"
LABEL software.version="0.2.0"
LABEL about.tags="Multiple Myeloma, ULP-WGS, ichorCNA"

RUN apt-get update && apt-get install -y \   
    libcurl4-gnutls-dev \
    libssl-dev \
    libssh2-1-dev \
    procps 

RUN R -e "install.packages(c('BiocManager', 'optparse'))"
RUN R -e "BiocManager::install(c('HMMcopy', 'GenomicRanges', 'GenomeInfoDb'))"

ENTRYPOINT ["Rscript"]
