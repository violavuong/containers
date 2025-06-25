FROM ubuntu:24.10
ARG buildDate

LABEL about.image="Ubuntu base image to run VarScan2 for SNV identification"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:24.10"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="VarScan2 image to detect SNV"
LABEL software.version="2.3.9"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \
	openjdk-17-jdk \
	wget
	
RUN mkdir varscan/ && \
	cd varscan/ && \
	wget https://sourceforge.net/projects/varscan/files/VarScan.v2.3.9.jar

ENTRYPOINT [ "java", "-jar", "/varscan/VarScan.v2.3.9.jar" ]
