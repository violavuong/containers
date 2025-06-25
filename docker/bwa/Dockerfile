FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Ubuntu base image to run bwa for sequence alignment"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="bwa to align fastq reads to the reference genome"
LABEL software.version="0.7.18"
LABEL about.tags="Multiple Myeloma, ULP-WGS, bwa"

RUN apt-get update && apt-get install -y \
    bwa \
    zlib1g-dev 

ENTRYPOINT ["bwa"]