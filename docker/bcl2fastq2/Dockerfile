FROM fedora:42
ARG buildDate

ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /tmp/

LABEL about.image="Fedora base image to run bcl2fastq2 for bcl-fastq conversion"
LABEL image.authors="Viola Vuong"
LABEL image.authors.email="violameixian.vuong2@unibo.it"
LABEL image.maintainer="Viola Vuong"
LABEL image.maintainer.email="violameixian.vuong2@unibo.it"
LABEL base.image="fedora:42"
LABEL image.version="1.0.0"
LABEL buildDate=${buildDate}
LABEL software="bcl2fastq2 for bcl-fastq conversion"
LABEL software.version="2.20"
LABEL about.tags="Multiple Myeloma, ULP-WGS, bcl2fastq2"

RUN sudo dnf -y install rpm-build unzip && dnf clean all
RUN unzip /tmp/bcl2fastq2-v2-20-0-linux-x86-64.zip && \
	sudo dnf -y install bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm
RUN chmod 777 /usr/local/bin/bcl2fastq

ENTRYPOINT [ "/usr/local/bin/bcl2fastq" ]

