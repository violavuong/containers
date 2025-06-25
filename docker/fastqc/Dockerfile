FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Linux Ubuntu base image to run vcflib to filter VCF files"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="fastqc image to manipulate fastq.gz files"
LABEL software.version="0.11.5"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \
	openjdk-23-jdk \
	perl \
	unzip

ENV FASTQC_LOC http://www.bioinformatics.babraham.ac.uk/projects/fastqc

WORKDIR /usr/local
ADD ${FASTQC_LOC}/fastqc_v0.11.5.zip /tmp/
RUN unzip /tmp/fastqc_v0.11.5.zip && sed -i 's/Xmx250m/Xmx2048m/' FastQC/fastqc && chmod 755 FastQC/fastqc
ENV PATH /usr/local/FastQC/:$PATH

RUN mkdir /data
WORKDIR /data

ADD ./wrapper.sh /usr/local/bin/
ENTRYPOINT ["bash", "/usr/local/bin/wrapper.sh"]
