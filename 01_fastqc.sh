#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=24,walltime=8:00:00
#PBS -A mcgaughs
#PBS -q mesabi
#PBS -N generating_fastqc_reports
#PBS -e /home/mcgaughs/jtmiller/process_all_popgen/error_out/$PBS_ARRAYID.$PBS_JOBID.err
#PBS -o /home/mcgaughs/jtmiller/process_all_popgen/error_out/$PBS_ARRAYID.$PBS_JOBID.out


#This takes every fastq file in the DIR and runs fastqc on it.
#	The home folder
export HOME=/home/mcgaughs/jtmiller/process_all_popgen
# The fastq
export DIR=/home/mcgaughs/shared/Datasets/{raw fastq)
# Load parallel
module load parallel

#	Define a bash function for doing the processing
#	since it is the same for each file, and we don't want to keep repeating

parfastqc() {
  #samples read from the glob below

 SAMP=$(grep "${1}" $HOME/metadata/fastq_sampleID.txt | awk '{print $2}')

 #run fastqc
 fastqc ${SAMP} --outdir=${HOME}/raw_fastqc_reports -q

}
#   Export function so we can call it with parallel
export -f parfastqc
#	cd into the reads directory

find $DIR -name *fq.gz | parallel parfastqc
