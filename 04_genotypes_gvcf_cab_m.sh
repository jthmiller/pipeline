#!/bin/bash -l
#PBS -A mcgaughs
#PBS -N genotypes_gvcfs_cab
#PBS -o /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/genotypes/gvcfs/er_out/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e /home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/genotypes/gvcfs/er_out/$PBS_JOBNAME.$PBS_JOBID.err

## Conda project environment
export PATH="$HOME/miniconda/bin:$PATH"
source activate cavefish_popgen

## Global environmental variables
## include option for java to designate the scratch location and the amount of
. /home/mcgaughs/jtmiller/popgen/cavefish_outliers/CODE/settings_gatk.sh $VER

# Set paths
## Location of output for large files  (ave. 23GB each)
outdir=/home/mcgaughs/shared/Datasets/outlier_analysis/gvcf_gz
## RAM drive. Load in compressed bam file for reading (ave. 12GB each)
##indir=/home/mcgaughs/jtmiller/popgen/cavefish_outliers/data/alignments
indir=/home/mcgaughs/shared/Datasets/outlier_analysis/bams
## Directory containing the fastqs that were aligned
fqdir=/home/mcgaughs/aherman/RIS_work/caballo_renamed_raw_fastqs/adapt_qual_trim
## ref director
refdir=/home/mcgaughs/jtmiller/amex_genomes

SCRDIR="/scratch.global/jtmiller"
TMPDIR="/scratch.local"

## ex: S1_CCGTCC_L23073_1_H5K73BCX2_surface.bam
### A will be L230 MOLNG1_S1
### B will be L289 MOLNG2_S1

samp=$(find $indir -name '*_L230[1-9]*_surface.bam' | xargs -n 1 basename | \
	awk -F'_' 'BEGIN{OFS="_"}{print $1,$2,$3}' | sort | uniq | sed -n "${PBS_ARRAYID}p")

echo "working on sample $samp"

bams=$(find $indir -name "${samp}*${VER}.bam")


#echo "copying bams to $TMPDIR"
#find $indir -name "$samp"_*_"$VER".bam -exec cp -v {} $TMPDIR/ \;
#ls -ltrh $TMPDIR
##samtools merge "$TMPDIR/temp_$samp.bam" "$TMPDIR"/"$samp"_*_"$VER".bam
###basebam=$(ls $TMPDIR/${samp}_*_${VER}.bam | xargs -n 1 basename)

### merge bams into RAMdisk #####
echo "merge bams in $TMPDIR"

samtools merge "$TMPDIR/${PBS_SERVER}_$samp.bam" $bams

du -h "$TMPDIR/${PBS_SERVER}_$samp.bam"

echo 'merged'

date

####ls "$TMPDIR"/"$samp"_*_"$VER".bam | xargs rm

### Index in memory

echo 'indexing...'

samtools index "$TMPDIR/${PBS_SERVER}_$samp.bam"

ls -lth "$TMPDIR/${PBS_SERVER}_$samp.bam"*

## Name appropriately
## ex: S1_CCGTCC_L23073_1_H5K73BCX2_surface.bam
### A will be L230 MOLNG1_S1
### B will be L289 MOLNG2_S1

RND=$(echo "$bams" | head -n1 | xargs -n 1 basename | grep -o _L2[0-9]*_ | sed 's/_//g')
RNDL=$(echo "$RND" | cut -c1-4)

if [ "$RNDL" == "L230" ]
then
	sampnew=$(echo $samp | sed 's/_Rnd2//g' |  awk -F'_' '{print "MOLNG1_"$1}')
else
	sampnew=$(echo $samp | sed 's/_Rnd2//g' |  awk -F'_' '{print "MOLNG2_"$1}')
fi

echo "this sample is $sampnew"

## Generate large .gvcf files from a single sample bam file (1 of 10 12GB files)
## GATK
module load gatk/3.7.0
which gatk

###
java "${HEAP}" -jar /panfs/roc/msisoft/gatk/3.7.0/GenomeAnalysisTK.jar\
	-T HaplotypeCaller\
	-R "${ref}"\
	-I "${TMPDIR}"/"${PBS_SERVER}"_"${samp}".bam\
	-o "${SCRDIR}"/"${PBS_SERVER}"_"${sampnew}".g.vcf\
	-nct "${NCT}"\
	--genotyping_mode DISCOVERY\
	--heterozygosity 0.005\
	--emitRefConfidence GVCF\
	-variant_index_type LINEAR\
	-variant_index_parameter 128000

wait

bgzip -f -c "${SCRDIR}/${PBS_SERVER}"_"${sampnew}".g.vcf > "${outdir}/${sampnew}".g.vcf.gz
tabix -p vcf "${outdir}/${sampnew}".g.vcf.gz

wait

## Cleanup
##rm $SCRDIR/"${PBS_SERVER}"_${sampnew}.g.vcf.gz
rm "$TMPDIR/"${PBS_SERVER}"_${samp}.bam"
rm "$TMPDIR/"${PBS_SERVER}"_${samp}.bam.bai"

echo -n "Done: "
date
