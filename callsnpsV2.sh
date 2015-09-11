#!/bin/bash

#$ -cwd
#$ -S /bin/bash
#$ -N snp4
#$ -o snpout
#$ -e snperr
# #$ -l mem_free=10G
#$ -V
#$ -p -10
#$ -t 1-4921:1
# #$ -t 4911-4921:1

i=$(expr $SGE_TASK_ID - 1)

PATH=$PATH:/raid1/home/bpp/knausb/bin/samtools-0.1.18/bcftools/

#REF="/raid1/labs/Grunwald_Lab/web_data/pinf_broad/mt/phytophthora_infestans_mito_haplotype_iia_1_supercontigs.fasta"
REF="/nfs0/labs/Grunwald_Lab/home/everhars/pram/ramorum1.fasta"

echo -n "Running on: "
hostname
echo "SGE job id: $JOB_ID"
date

# BAMS must be sorted.
# samtools mpileup
# -u uncompressed bcf
# -g compute genotype likelihoods
# -f faidx fasta

# bcftools view
# -b output bcf [default vcf]
# -c Call variant using Bayesian inference
# -v output variant sites only
# -g Call per-sample genotypes at variant sites

#samtools mpileup -P ILLUMINA -ugf $REF ../*/*sorted.bam | bcftools view -bcvg - > var.raw.bcf
#samtools mpileup -P ILLUMINA -ugf $REF ../*/*sorted.bam > pinf.pileup.bcf
# | bcftools view -bcvg - > var.raw.bcf

#samtools mpileup -P ILLUMINA -ugf $REF -b bams.txt > pinf.pileup.bcf
#
#/home/bpp/knausb/bin/samtools-0.1.18/bcftools/bcftools view -bcvg pinf.pileup.bcf > pinf2.raw.bcf
#
#
# samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf

#CMD="samtools mpileup -P ILLUMINA -ugf $REF -b bams.txt -r Supercontig_1."$SGE_TASK_ID" > ./bcfs/pinf_1.$SGE_TASK_ID.pileup.bcf"
#CMD="samtools mpileup -ugf $REF -b bams.txt -r Supercontig_1."$SGE_TASK_ID" | bcftools view -bvcg - > ./bcfs/pinf_1."$SGE_TASK_ID".raw.bcf"
#CMD="samtools mpileup -ugf $REF -b bams.txt -r Supercontig_1."$SGE_TASK_ID" - > pinf_1."$SGE_TASK_ID".pileup.bcf"
# | bcftools view -bvcg - > ./bcfs/pinf_1."$SGE_TASK_ID".raw.bcf"
#CMD="samtools mpileup -ugf $REF -b bams.txt -r Supercontig_1.$SGE_TASK_ID"

samtools mpileup -ugf $REF -b bams.txt -r Supercontig_1.$SGE_TASK_ID | bcftools view -bvcg - > ./bcfs/pinf_1.$SGE_TASK_ID.raw.bcf


# | bcftools view -bvcg - > ./bcfs/pinf_1."$SGE_TASK_ID".raw.bcf"
#echo $CMD
# $CMD
#| bcftools view -bcvg > bcfs/pinf_1.$SGE_TASK_ID.raw.bcf"


date

# vcfutils.pl varFilter
# -d INT minimum read depth [2]
# -D INT maximum read depth [10000000]

#bcftools view var.raw.bcf | vcfutils.pl varFilter -D 2000 > var.flt.vcf
#bcftools view pinf2.raw.bcf | vcfutils.pl varFilter -D 2000 > var.flt.vcf

#date

# EOF.
