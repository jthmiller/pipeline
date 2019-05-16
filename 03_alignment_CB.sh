#!/bin/bash -l
#PBS -A mcgaughs
#PBS -e /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments/er_out/$PBS_ARRAYID.$PBS_JOBID.err
#PBS -o /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments/er_out/$PBS_ARRAYID.$PBS_JOBID.out

## Project environment
export PATH="$HOME/miniconda/bin:$PATH"
source activate cavefish_popgen

## directories
### indir=/home/mcgaughs/aherman/RIS_work/caballo_renamed_raw_fastqs
### indir=/home/mcgaughs/shared/Datasets/outlier_analysis/aligned_fq
### indir=/home/mcgaughs/shared/Datasets/outlier_analysis/aligned_fq

indir=/home/mcgaughs/shared/Datasets/RAW_NGS/Caballo_Moro_trimmed_fq
outdir=/home/mcgaughs/shared/Datasets/bams/v1_cavefish_Caballo_Moro
#### outdir=/home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments
## outdir=/home/mcgaughs/shared/Datasets/outlier_analysis/bams ### change to this
## fqdir=/panfs/roc/groups/14/mcgaughs/grossj
refdir=/home/mcgaughs/jtmiller/amex_genomes

## Global environmental variables
. /home/mcgaughs/jtmiller/popgen/cavefish_outliers/CODE/settings_gatk.sh $VER
# define $VER if running this script alone in the qsub command for submission
indir=/home/mcgaughs/shared/Datasets/Reads_ready_to_align/Caballo_Moro

echo $PBS_ARRAYID

## Finds all the forward fq in DIR _adtrim_trim_pair_R1
## fq1=$(ls "$shar"/*_adtrim_trim_pair_R1.fastq | sed -n "${PBS_ARRAYID}p") ## just for missed alignments

fq1=$(ls "$indir"/*_L2*_adtrim_trim_pair_R1.fastq.gz | sed -n "${PBS_ARRAYID}p")
##### fq1=$(find $indir -name *_adtrim_trim_pair_R1.fastq.gz | sed -n "${PBS_ARRAYID}p")
fq2=$(echo $fq1 | sed 's/_R1/_R2/')
#### root sample name

echo $fq1
echo $fq2

###names
root=$(basename $fq1 | sed 's/_adtrim_trim_pair_R1\.fastq\.gz//')
echo $root

## Save fq header
FQID=$(zcat "$fq1" | head -n1)
## Read group of each fastq
####SAMP=$(echo $root | awk -F"_" '{print $1}')

SAMP=$(echo $root | awk -F"_" '{print $1"_"$2}')
INST=$(echo "$FQID" | awk -F":" 'BEGIN{OFS="\t"}{print $1}') ## Insturment ID
RUN=$(echo "$FQID" | awk -F":" 'BEGIN{OFS="\t"}{print $2}') ## internal run number
FCID=$(echo "$FQID" | awk -F":" 'BEGIN{OFS="\t"}{print $3}') ## Flowcell ID, ex: H5K73BCX2
LANE=$(echo "$FQID" |awk -F":" 'BEGIN{OFS="\t"}{print $4}') ## Lane ID 1 or 2

echo "SAMP" "INST" "RUN" "FCID" "LANE"
echo "$SAMP" "$INST" "$RUN" "$FCID" "$LANE"

## Set RG info
##### SBC=$(echo "$root" | awk -F"_" 'BEGIN{OFS="\t"}{print $2}') ## Sample barcode

SBC=$(echo "$root" | awk -F"_" 'BEGIN{OFS="\t"}{print $2}') ## Sample barcode
ID=$(echo "$FCID.$LANE") ## {FLOWCELL_BARCODE}.{LANE}
###LB=$(echo "$root" | awk -F"_" 'BEGIN{OFS="\t"}{print $3}') ## Run 1 L230XX is MOLNG-1690, Run 2 L289XX is MOLNG-2003
LB=$(echo "$root" | awk -F"_" 'BEGIN{OFS="\t"}{print $4}') ## Run 1 L230XX is MOLNG-1690, Run 2 L289XX is MOLNG-2003
PL="Illumina"
PU=$(echo "$FCID.$LANE.$SBC") ##  {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}
RG="@RG\\tID:$ID\\tPL:$PL\\tPU:$PU\\tLB:$LB\\tSM:$SAMP"

echo "READ GROUP"
echo $RG

## scratch and out file names
outfile=${outdir}/${root}_$VER.bam
echo $outfile

## going to work
###bwa mem -t "$NCT" -k 12 -M -R $RG $ref $fq1 $fq2 | \ ### Alignment
###samblaster -M | \ #### dups/split reads with -M for compatibility with picard
###samtools view -hbu - | samtools sort -m $MEG -@ $NCT -o - -T $SCRDIR/$root.tmp > $outfile

bwa mem -t "${NCT}" -k 12 -M -R $RG $ref $fq1 $fq2 | ### Alignment
samblaster -M |  #### dups/split reads with -M for compatibility with picard
samtools view -hbu - > $SCRDIR/$root.unsorted.bam
samtools sort -m $MEG -@ $NCT -o - -T $SCRDIR/$root.tmp $SCRDIR/$root.unsorted.bam > $outfile

echo "Done"
