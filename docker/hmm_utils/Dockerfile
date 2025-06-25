FROM ubuntu:24.10
ARG buildDate

ADD /hmmcopy_utils/ /hmmcopy_utils/

LABEL about.image="Ubuntu base image to run HMMcopy to create wig files from bam files"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="HMMcopy utils to convert bam files into wig files"
LABEL software.version="0.99.0"
LABEL about.tags="Multiple Myeloma, ULP-WGS, ichorCNA"

RUN apt-get update && apt-get install -y \
    build-essential
