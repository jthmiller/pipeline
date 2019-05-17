# pipeline

### To run an interactive session on Itasca
qsub -I -q lab -lwalltime=8:00:00,nodes=1:ppn=6,mem=24gb



## Genotype Calling
GATK4 Haplotypecaller with '--emitRefConfidence GVCF' to get reference allele depth and confidence. Results in individual gvcfs.
GATK4 GenomicsDB makes a database for calling SNPs. Takes the longest, and requires all samples present. Must be redone if sample added.
GATK4 GenotypeGVCFs uses database for calling SNPs on individual GVCFs.
