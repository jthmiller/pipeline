#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=24,walltime=2:00:00
#PBS -A mcgaughs
#PBS -m abe
#PBS -j oe
#PBS -M jtmiller@umn.edu
#PBS -q mesabi
#PBS -N generating_fastqc_reports

#	The home folder
export HOME=/home/mcgaughs/jtmiller/popgen/cavefish_RIS/
# The fastq
export DIR=/panfs/roc/groups/14/mcgaughs/grossj/

module load parallel/20160822
module load fastqc/0.11.7

#	Define a bash function for doing the processing
#	since it is the same for each file, and we don't want to keep repeating

parfastqc() {
	#samples read from the glob below
	SAMP=${1}

	#run fastqc
	fastqc ${SAMP} -o ${HOME}raw_fastqc_reports -q
}
#   Export function so we can call it with parallel
export -f parfastqc
#	cd into the reads directory

find $DIR -name *fq.gz | parallel parfastqc


###parallel parfastqc ::: "find $DIR -name *fq.gz"

###cd ${HOME}caballo_renamed_raw_fastqs
###parallel parfastqc ::: *.fastq.gz
