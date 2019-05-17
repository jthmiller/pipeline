#!/bin/bash -l
#PBS -A mcgaughs
#PBS -e /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments/er_out/$PBS_ARRAYID.$PBS_JOBID.err
#PBS -o /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments/er_out/$PBS_ARRAYID.$PBS_JOBID.out

## Project environment
export PATH="$HOME/miniconda/bin:$PATH"
source activate cavefish_popgen

export IN_DIR
export OUT_DIR
export HOME

# define $VER if running this script alone in the qsub command for submission
. /home/mcgaughs/jtmiller/popgen/cavefish_outliers/CODE/settings_gatk.sh $VER
export $ref

echo $PBS_ARRAYID

smpl=$(sed -n "${PBS_ARRAYID}p" $SAMPMAP | awk '{print $2}')

## Paired reads
fq1=$(ls "${IN_DIR}"/*_trim_pair_R1.fastq.gz | grep $smpl)
##### fq1=$(find $indir -name *_adtrim_trim_pair_R1.fastq.gz | sed -n "${PBS_ARRAYID}p")
fq2=$(echo $fq1 | sed 's/_R1/_R2/')

## Unpaired reads
fq1up=$(echo $fq1 | sed 's/_trim_pair_/_trim_unpair_/')
fq2up=$(echo $fq2 | sed 's/_trim_pair_/_trim_unpair_/')

#### root sample name
## the fq2 calls the fq1 and simply replaces "R1" with "R2"
echo $fq1
echo $fq2

###names
root=$(basename $fq1 | sed 's/_adtrim_trim_pair_R1\.fastq\.gz//')
echo $root
#we also need to name these as population_sample#
## Save fq header
FQID=$(zcat "$fq1" | head -n1)
## Read group of each fastq
####SAMP=$(echo $root | awk -F"_" '{print $1}')

## still renaming below
SAMP=$(echo $root | awk -F"_" '{print $1"_"$2}')
INST=$(echo "$FQID" | awk -F":" 'BEGIN{OFS="\t"}{print $1}') ## Insturment ID
RUN=$(echo "$FQID" | awk -F":" 'BEGIN{OFS="\t"}{print $2}') ## internal run number
FCID=$(echo "$FQID" | awk -F":" 'BEGIN{OFS="\t"}{print $3}') ## Flowcell ID, ex: H5K73BCX2
LANE=$(echo "$FQID" |awk -F":" 'BEGIN{OFS="\t"}{print $4}') ## Lane ID 1 or 2

echo "SAMP" "INST" "RUN" "FCID" "LANE"
echo "$SAMP" "$INST" "$RUN" "$FCID" "$LANE"


#?should be zcat line for parsing? JEFF TODO?
## Set RG info
##### SBC=$(echo "$root" | awk -F"_" 'BEGIN{OFS="\t"}{print $2}') ## Sample barcode
## GET FASTQ HEADER
FQHEADER=$(zgrep -m 1 '^@' "$fq1")

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
#
out_pair=${outdir}/${root}_${VER}_paired_sorted.bam
out_unpaired=${outdir}/${root}_${VER}_sorted.bam

echo $out_pair
echo $out_unpair

#### dups/split reads with -M for compatibility with picard

### PAIRED READS ###############################################################
bwa mem -t "${NCT}" -k 12 -M -R $RG $ref $fq1 $fq2 \
| samblaster -M \
| samtools view -hbu - \
| samtools sort -m $MEG -@ $NCT -o - -T $SCRDIR/$root.tmp > $out_pair
################################################################################

### UNPAIRED ###################################################################
cat $fq1up $fq2up ${SCRDIR}/${root}_unpaired.fq

bwa mem -t "${NCT}" -k 12 -M -R $RG $ref ${SCRDIR}/${root}_unpaired.fq \
| samblaster -M \
| samtools view -hbu - \
| samtools sort -m $MEG -@ $NCT -o - -T $SCRDIR/$root.tmp > $out_unpair
################################################################################

echo "Done"
