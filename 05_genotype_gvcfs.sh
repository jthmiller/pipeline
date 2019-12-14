#!/bin/bash -l
#PBS -A mcgaughs
#PBS -N genotypes_vcf_intervals
#PBS -o /scratch.global/jtmiller/er_out/$PBS_JOBNAME.$PBS_JOBID.out
#PBS -e /scratch.global/jtmiller/er_out/$PBS_JOBNAME.$PBS_JOBID.err

## Conda project environment
export PATH="$HOME/miniconda/bin:$PATH"
source activate cavefish_popgen
SHARE=/home/mcgaughs/shared/Datasets/Miller_outlier_analysis

### Set paths to shared data
indir=$SHARE/gvcf_gz
outdir=$SHARE/interval_gvcf
metadir=$SHARE/metadata
refdir=$HOME/amex_genomes

## Global environmental variables
. /home/mcgaughs/jtmiller/popgen/cavefish_outliers/CODE/settings_gatk.sh $VER

## interval is half a megabase
chr=$(sed -n "${PBS_ARRAYID}p" "${ref%\.fna}".windows.500KB.tsv | awk '{print $1":"$2"-"$3}')
###~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.tsv
## outcheck
echo $PBS_O_QUEUE
echo $chr

## SWITCHING TO GATK 4
## Databases done
/home/mcgaughs/jtmiller/software/gatk-4.1.0.0/gatk --java-options "${HEAP}" GenomicsDBImport \
 --tmp-dir="${TMPDIR}" \
 --sample-name-map "${metadir}"/cohort.sample_map \
 --genomicsdb-workspace-path "${SCRDIR}"/DBs/scaffolds/"${chr}" \
 --batch-size 60 \
 -L "${chr}" \
 --reader-threads "${NCT}" \
 --overwrite-existing-genomicsdb-workspace

### Echo
ls "${SCRDIR}"/DBs/scaffolds/"${chr}"

du -sh "${SCRDIR}"/DBs/scaffolds/"${chr}"

echo "done with database"

/home/mcgaughs/jtmiller/software/gatk-4.1.0.0/gatk --java-options "${HEAP}" GenotypeGVCFs \
 -R "${ref}" \
 --tmp-dir="${TMPDIR}" \
 -all-sites true \
 -V gendb://"${SCRDIR}"/DBs/scaffolds/"${chr}" \
 -G StandardAnnotation -new-qual \
 -O "${SCRDIR}/scaff_gvcf/${chr}.vcf" \
 --heterozygosity 0.005 \
 -L "${chr}"

### bgzip files bc they are huge
bgzip -c "${SCRDIR}/scaff_gvcf/${chr}.vcf" > "${outdir}/${chr}.vcf.gz"
tabix -p vcf "${outdir}/${chr}.vcf.gz"
bcftools index "${outdir}/${chr}.vcf.gz"

echo -n "Done: "
date
