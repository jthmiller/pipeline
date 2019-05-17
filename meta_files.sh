
##### Get v2_surface refernce ##################################################
################################################################################
curl ftp://ftp.ensembl.org/pub/release-96/fasta/astyanax_mexicanus/dna/Astyanax_mexicanus.Astyanax_mexicanus-2.0.dna.toplevel.fa.gz \
 -o "${refdir}"/Astyanax_mexicanus.Astyanax_mexicanus-2.0.dna.toplevel.fa.gz
gunzip /Astyanax_mexicanus.Astyanax_mexicanus-2.0.dna.toplevel.fa.gz
################################################################################


#### SAMPLE MAP FOR GATK #######################################################
METADIR=/home/mcgaughs/shared/Datasets/outlier_analysis/metadata
gvcfdir=/home/mcgaughs/shared/Datasets/outlier_analysis/gvcf_gz
find $gvcfdir -maxdepth 1 -name '*vcf.gz' > ${METADIR}/cohort.sample_map_b
cat "${METADIR}"/cohort.sample_map_b | xargs -n 1 basename | sed 's/.g.vcf.gz//g' > "${METADIR}"/cohort.sample_map_a
paste -d '\t' "${METADIR}"/cohort.sample_map_a "${METADIR}"/cohort.sample_map_b > "${METADIR}"/cohort.sample_map
#### SAMPLE MAP FOR GATK #######################################################
################################################################################

#### Make intervals for GATK ###################################################
#### 500kb windows #############################################################
bedtools makewindows -b ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.bed \
-w 500000 > ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.bed

awk '{print $1,$2,$3}' ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.bed > ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.bed.tsv

mv ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.bed.tsv ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.tsv.bed

echo -e "contig\t" "start\t" "stop" > ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.tsv
awk '{print $1,($2+1),$3}' ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.tsv.bed >> ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.500KB.tsv

#### 200kb windows #############################################################
bedtools makewindows -b ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.bed \
-w 200000 > ~/amex_genomes/GCF_000372685.2_Astyanax_mexicanus-2.0_genomic.windows.200KB.bed
### Make intervals for GATK ###################################################
################################################################################

### make list of non-empty interval vcfs #######################################
################################################################################
 inters=$(awk '{print $1":"$2"-"$3}' "${ref%\.fna}.windows.200KB.bed")
 for i in "$iners" ; do find $snpsindir -type f -size +1c -name "$i"*.vcf.gz ; done > $metadir/interval_vcf_200K.list
######################################################
### Index intervals ##################################
 while read l; do bcftools index "${l}" ; done < $metadir/interval_vcf_200K.list
### make list of non-empty interval vcfs #######################################
################################################################################

## Make a BCF of all intervals ######################################################
#####################################################################################
bcftools concat --threads $NCT -a --file-list $metadir/interval_vcf_200K.list -O b -o "${outdir}"/"${VER}"_snps.bcf
bcftools index "${outdir}"/"${VER}"_snps.bcf
metadir=/home/mcgaughs/shared/Datasets/outlier_analysis/metadata
bcftools concat -a --file-list $metadir/interval_vcf_200K.list | head
bcftools sort -O b "${outdir}"/"${VER}"_snps.bcf -o try.bcf
chmod 555 "${outdir}"/"${VER}"_snps.bcf
## Make a BCF of all intervals ######################################################
################################################################################

##### Generate Depth Stats for hard filter (take out extremely low sites)############
#####################################################################################
bcftools query -e'MIN(INFO/DP)<10' -r 'NC_035899.1' -f'[%CHROM:%POS\t%INFO/QD\t%INFO/DP\n]' "${outdir}"/"${VER}"_snps.bcf | gzip > $metadir/"${VER}"_snps.depth
qsub -t 1 -v VER="surface", NCT="12" -l mem=64gb -l nodes=1:ppn=12 -l walltime=24:00:00 -q cavefish ${CO}/CODE/07_filter_variants.sh
#####################################################################################
### Check out depth of samples and make a hsitogram
bcftools stats --depth 100,6000,50 -e 'INFO/DP<50' "${outdir}"/"${VER}"_snps.bcf > "${metadir}"/all.depth.stats
## falls off fast at 1-25 (x50) = ~1500
#####################################################################################

### PHASE GENOTYPES##################################################################
