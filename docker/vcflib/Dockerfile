FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Ubuntu base image to run vcflib for filtering VCF files"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="vcflib image to filter VCF files"
LABEL software.version="1.0.10"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \
	libvcflib-tools \
	libvcflib-dev

