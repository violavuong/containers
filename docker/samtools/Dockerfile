FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Ubuntu base image to run samtools for sorting and indexing aligned fastq files"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="samtools to process fastq files"
LABEL software.version="1.20"
LABEL about.tags="Multiple Myeloma, ULP-WGS, samtools"

RUN apt-get update && apt-get install -y \
    gcc \
    samtools

ENTRYPOINT ["samtools"]
