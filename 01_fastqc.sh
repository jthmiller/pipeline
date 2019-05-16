#!/bin/bash
#PBS -l mem=62gb,nodes=1:ppn=24,walltime=8:00:00
#PBS -A mcgaughs
#PBS -q mesabi
#PBS -N generating_fastqc_reports
#PBS -e /home/mcgaughs/jtmiller/process_all_popgen/error_out/$PBS_ARRAYID.$PBS_JOBID.err
#PBS -o /home/mcgaughs/jtmiller/process_all_popgen/error_out/$PBS_ARRAYID.$PBS_JOBID.out

export IN_DIR
export OUT_DIR
export HOME

#This takes every fastq file in the DIR and runs fastqc on it.
#	The home folder

# Load parallel
module load parallel

#	Define a bash function for doing the processing
#	since it is the same for each file, and we don't want to keep repeating

parfastqc() {
  #samples read from the glob below
 #run fastqc
 fastqc ${1} --outdir=${HOME}/raw_fastqc_reports -q

}
#   Export function so we can call it with parallel
export -f parfastqc
#	cd into the reads directory

find $IN_DIR -type f -name *fq.gz | parallel parfastqc
