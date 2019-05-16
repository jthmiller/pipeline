#!/bin/bash
################################################################################
## ENVIRONMENT #################################################################
## Create a conda environment containing software ##############################
## versions used in the Adam et al pipeline ####################################

CONDA=/home/mcgaughs/jtmiller/bin/miniconda3/bin/conda

if [ cavefish_popgen env not avail ]
then
	if [ CONDA is null]
	echo "install conda"
	die
fi
## Use own version of conda
alias conda=$CONDA
## set up channels for bioconda programs
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

### Conda env with versioned software
conda create -n cavefish_popgen samblaster fastqc bwa==0.7.4 samtools picard==2.3.0 vcftools bedtools plink tabix gatk4=4.1.0.0 bcftools
source activate cavefish_popgen
fi
################################################################################
################################################################################
### Put output directory here
HOME=/home/mcgaughs/jtmiller/process_all_popgen
SHARED=

mkdir $HOME/error_out
mkdir $HOME/metadata
mkdir $HOME/raw_fastqc_reports

################################################################################
## V2_Surface Genome ###########################################################

qsub \
 -l mem=120gb \
 -l nodes=1:ppn=24 \
 -l walltime=8:00:00 \
 -q cavefish \
 -v HOME=/home/mcgaughs/jtmiller/process_all_popgen \
 -v OUT_DIR=$HOME/raw_fastqc_reports \
 -v IN_DIR=$SHARED/ (RAW FASTQ DIR) \
 ${HOME}/pipeline/01_fastqc.sh

### Requires sample name map. Tab delimited two column file with "*.fastq   {sample_name}"
qsub \
 -t 1-8 \
 -l mem=21gb \
 -l nodes=1:ppn=24 \
 -l walltime=8:00:00 \
 -q cavefish \
 -v SAMPMAP=$HOME/metadata/sample_map.txt \
 -v HOME=/home/mcgaughs/jtmiller/popgen/ \
 -v OUT_DIR=$SHARED/(shared trimmed fq folder) \
 -v IN_DIR=$SHARED/(shared raw fq folder) \
 ${HOME}/pipeline/02_quality_trim_CB_B.sh

qsub \
 -t 1-23 \
 -v VER="surface" \
 -N alignment \
 -l mem=62gb \
 -l nodes=1:ssd:ppn=12 \
 -l walltime=48:00:00 \
 -q cavefish
 -v HOME=/home/mcgaughs/jtmiller/popgen/ \
 -v IN_DIR= (shared trimmed fq folder) \
 -v OUT_DIR= (shared bam dir) \
 -v VER=surface \
 ${HOME}/pipeline/03_alignment.sh)


indir=/home/mcgaughs/shared/Datasets/RAW_NGS/Caballo_Moro_trimmed_fq
outdir=/home/mcgaughs/shared/Datasets/bams/v1_cavefish_Caballo_Moro
#### outdir=/home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments
## outdir=/home/mcgaughs/shared/Datasets/outlier_analysis/bams ### change to this
## fqdir=/panfs/roc/groups/14/mcgaughs/grossj
refdir=/home/mcgaughs/jtmiller/amex_genomes






### Generate Gvcfs
qsub \
 -t 1-13 \
 -v VER="surface",SCRDIR=/scratch.ssd,TMPDIR=/dev/shm,_JAVA_OPTIONS="-Djava.io.tmpdir=/dev/shm -Xmx15g",NCT=12 \
 -l mem=31gb \
 -l nodes=1:ssd:ppn=12 \
 -l walltime=14:00:00 \
 -q mesabi \
${CO}/CODE/04_genotypes_gvcf_jineo_escon.sh)

#### Make a sample map for GATK (see meta_files.sh #############################

#### Genotype GVCFs with GenomicDB #############################################
## Total intervals: 4668, 1 tourqe job per interval ############################

qsub -t 1-1750 \
 -v VER="surface",NCT="1",HEAP='-Xmx6G',SCRDIR='/scratch.global/jtmiller',TMPDIR='/scratch.local',_JAVA_OPTIONS='-Xmx6G' \
 -l mem=8gb -l nodes=1:ppn=2 -l walltime=24:00:00 \
 -q mesabi ${CO}/CODE/06_genotypes_vcf.sh

################################################################################
################################################################################


### CAVE GENOME ALIGNMENT ######################################################
################################################################################
qsub \
 -t 1-23 \
 -v VER="cave",NCT='12',MEG='28G',SCRDIR='/scratch.global/alignments' \
 -N CM_alignment \
 -l mem=32gb \
 -l nodes=1:ppn=12 \
 -l walltime=48:00:00 \
 -q mesabi \
 ${CO}/CODE/03_alignment_b.sh

qsub \
 -t 1-40 \
 -v VER="cave",NCT='12',MEG='28G',SCRDIR='/scratch.global/alignments' \
 -N CB_alignment \
 -l mem=32gb \
 -l nodes=1:ppn=12 \
 -l walltime=24:00:00 \
 -q mesabi \
 ${CO}/CODE/03_alignment_CB_b.sh
################################################################################
################################################################################
