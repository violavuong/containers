FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Ubuntu base image to run gatk3 on bam files"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="gatk 3.8.1"
LABEL software.version="3.8.1"
LABEL about.tags="Multiple Myeloma, ULP-WGS, gatk3"

RUN apt-get update && apt-get install -y \
    lbzip2 \
    openjdk-8-jdk \
    tar \
    wget
	
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
    tar -xvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
    mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk3.8-1-0/

ENTRYPOINT ["java", "-jar", "/gatk3.8-1-0/GenomeAnalysisTK.jar", "-T"]