FROM bzhanglab/copywriter:2.20.0
ARG buildDate

ADD /copywriteR/ /cwR_lib/

LABEL about.image="R base image to run CopywriteR for CNAs identification"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="copywriter:2.20.0"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="CopywriteR image to detect CNA"
LABEL software.version="2.20.6"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \
	libcurl4-gnutls-dev \
	libssl-dev \
	libssh2-1-dev \
	procps

RUN R -e 'install.packages("optparse")'

ENTRYPOINT [ "Rscript" ]
