FROM ubuntu:25.04
ARG buildDate

LABEL about.image="Ubuntu base image to run fastp"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="ubuntu:25.04"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="fastp for FASTQ pre-processing"
LABEL software.version="1.0.1"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"


RUN apt-get update && apt-get install -y \
    wget
	
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp
	
ENTRYPOINT ["./fastp"]
