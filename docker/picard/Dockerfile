FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Ubuntu base image to run picard for duplicates removal"
LABEL image.authors="Gaia Mazzocchetti, Viola Vuong"
LABEL image.authors.email="gaia.mazzocchetti2@unibo.it, violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="picard to process fastq files"
LABEL software.version="3.3.0"
LABEL about.tags="Multiple Myeloma, ULP-WGS, picard"

RUN apt-get update && apt-get install -y \
    openjdk-17-jdk \
    git 
	
RUN git clone https://github.com/broadinstitute/picard.git 
RUN cd /picard && ./gradlew shadowJar

ENTRYPOINT ["java", "-jar", "/picard/build/libs/picard.jar"]