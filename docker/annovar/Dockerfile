FROM perl:5.41.5
ARG buildDate

LABEL about.image="Perl base image to run Annovar for variants annotation"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="perl:5.41.5"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="Annovar to annotate variants"
LABEL software.version="20.06.08"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \
	tar \
	wget

RUN wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz && \
	tar -xvzf annovar.latest.tar.gz && \
	rm annovar.latest.tar.gz

ENV LC_CTYPE=en_US.UTF-8
ENV LC_ALL=en_US.UTF-8
