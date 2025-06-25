#!/bin/bash

set -e

print_usage(){
>&2 cat <<EOF
This docker image runs fastq on an input directory with .fastq or .fastq.gz files.
Output is put in the same directory consists of a <input ID>_fastq.zip file and a <input ID>_fastq.html file
which shows an graphical overview of results.

For details see http://www.bioinformatics.babraham.ac.uk/projects/fastqc
EOF
}

help=$1
if [ ! -z "$1" ]; then
    print_usage
    exit
fi

finish() {
    # Fix ownership of output files
    uid=$(stat -c '%u:%g' /data)
    chown -R $uid /data
}
trap finish EXIT

for file in /data/*.fastq
do
	fastqc $file
done

for file in /data/*.fastq.gz
do
	fastqc $file
done

exit
