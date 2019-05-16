#!/bin/bash
#PBS -l mem=22gb,nodes=1:ppn=8,walltime=24:00:00
#PBS -A mcgaughs
#PBS -N par_trim
#PBS -q cavefish
#PBS -o /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/genotypes/gvcfs/er_out/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/genotypes/gvcfs/er_out/$PBS_JOBNAME.$PBS_JOBID.err

export IN_DIR
export OUT_DIR
export HOME
export SAMPMAP

TRIM_DIR=/panfs/roc/groups/14/mcgaughs/smcgaugh/tools/Trimmomatic-0.30
 #change to /panfs/roc/groups/14/mcgaughs/shared/Datasets/Reads_ready_to_align

module load parallel

# bash function for applying trimmomatic with gnu parallel
#output name convention: population_sample

partrim() {
 #sample basename from the paired_basenames.txt file

 BASE=$(basename ${1})

 newname=$(grep "${BASE}" "${SAMPMAP}" | awk '{print $2}')

 pair=$(echo $1 | sed 's/R1/R2/')

 #trim low quality
 java -Xmx2g -jar ${TRIM_DIR}/trimmomatic-0.30.jar PE -phred33 \
  ${1} ${pair} \
  ${OUT_DIR}/${newname}_trim_pair_R1.fastq.gz ${OUT_DIR}/${newname}_trim_unpair_R1.fastq.gz \
  ${OUT_DIR}/${newname}_trim_pair_R2.fastq.gz ${OUT_DIR}/${newname}_trim_unpair_R2.fastq.gz \
  SLIDINGWINDOW:6:30 MINLEN:50
}
#   Export function so we can call it with parallel
export -f partrim

 parallel --joblog ${OUT_DIR}/trimmomatic_parallel_logfile.txt partrim

echo 'DONE'

### Both paired and orphaned reads are aligned in 03_alignment
