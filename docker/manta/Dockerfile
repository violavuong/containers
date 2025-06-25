FROM python:2.7.17
ARG buildDate

LABEL about.image="Python2 base image to run Manta for t-IgH identification"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="python:2.7.17"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="Manta image to detect t-IgH"
LABEL software.version="1.6.0"
LABEL about.tags="Multiple Myeloma, UMA panel, bioinformatics"

RUN apt-get update && apt-get install -y \
	build-essential \
	bzip2 \
	gcc \
	g++ \
	make \
	tar \
	wget \
	zlib1g-dev

RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.release_src.tar.bz2 && \
	tar -xjf manta-1.6.0.release_src.tar.bz2 && \
	mkdir build && cd build && \
	../manta-1.6.0.release_src/configure --jobs=4 --prefix=/build/src/python/ && \
	make -j4 install

ENTRYPOINT [ "python" ]
