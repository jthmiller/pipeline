#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -A mcgaughs
#PBS -m abe
#PBS -j oe
#PBS -M jtmiller@umn.edu
#PBS -q batch
#PBS -N par_trim
## 12 cores 128GB PBS -q mcgaugh
#PBS -o /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/genotypes/gvcfs/er_out/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/genotypes/gvcfs/er_out/$PBS_JOBNAME.$PBS_JOBID.err


HOME=/home/mcgaughs/jtmiller/popgen/
TRIM_DIR=/panfs/roc/groups/14/mcgaughs/smcgaugh/tools/Trimmomatic-0.30/
OUT_DIR=/home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/qual_trim_fq/
DIR=/panfs/roc/groups/14/mcgaughs/grossj/
# modules
module load parallel/20160822
module load java/jdk1.7.0_45

# list of samples put in outdir
find $DIR -name *fq.gz | awk '{print $9}' FS="\/|_" | sed 's/[1:2].fq.gz//g' | uniq > ${OUT_DIR}Esc_Jin_samp_List.txt
# bash function for applying trimmomatic with gnu parallel

partrim() {
	#sample basename from the paired_basenames.txt file
	BASE=${1}
	SAMP1=$(find $DIR -type f -name ${BASE}1.fq.gz)
	SAMP2=$(find $DIR -type f -name ${BASE}2.fq.gz)

	#trim low quality
	java -Xmx2g -jar ${TRIM_DIR}trimmomatic-0.30.jar PE -phred33\
		${SAMP1} ${SAMP2}\
		${OUT_DIR}${BASE}trim_pair_R1.fastq.gz ${OUT_DIR}${BASE}trim_unpair_R1.fastq.gz\
		${OUT_DIR}${BASE}trim_pair_R2.fastq.gz ${OUT_DIR}${BASE}trim_unpair_R2.fastq.gz\
		SLIDINGWINDOW:6:30 MINLEN:50
}
#   Export function so we can call it with parallel
export -f partrim

parallel --joblog ${OUT_DIR}trimmomatic_parallel_logfile.txt -a ${OUT_DIR}Esc_Jin_samp_List.txt partrim
