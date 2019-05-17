### Make a plink ID file from plinks output ####################################
################################################################################
awk -F'_| ' '{print $1,$1"_"$2,$1"_cave"}' surface_masked_indel3.test.fam > surface_masked_indel3.within
plink --vcf "${outdir}"/"${VER}"_masked_indel3.test.vcf.gz --out ${plinkdir}/"${VER}"_masked_indel3.test \
 --make-bed --allow-extra-chr --autosome-num 24 --allow-no-sex --double-id \
 --within ~/popgen/cavefish_outliers/notes/surface_masked_indel3.within --write-cluster

awk -F'_| ' '{print $1,$1"_"$2,$1"_cave"}' surface_masked_indel3.test.fam > surface_masked_indel3.within
awk -F' ' '{print "0",$2,$3}' ~/popgen/cavefish_outliers/notes/surface_masked_indel3.within > ~/popgen/cavefish_outliers/notes/surface_masked_indel3.within.contID
## Make a plink ID file from plinks output #####################################
################################################################################
