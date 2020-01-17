#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=01:25:00

# Output
#$ -j y

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
source /broad/software/scripts/useuse

# Use your dotkit
reuse Bcftools

##################
### Run script ###
##################

# Aim 3, compare allele fraction for WGS bam to ATAC bam for vcf variants from Archer 2018 

# $1 The sample to be analyzed
# $2 [optional] if the ATAC bam file is named something else, add it as a second argument.

SAMPLE=$1
ABAM=${2:-$1}
ATAC=/broad/medullo/ATACseq_fraenkel

mkdir -p $ATAC/allele_fraction

PREFIX=$ATAC/allele_fraction/${SAMPLE}_AF
echo "Doing mpileup on sample $SAMPLE..."

START=$SECONDS
bcftools mpileup -f $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	-T /broad/medullo/VCF_proteomics_cohort/${SAMPLE}_pair.pass.vcf\
	-a AD\
	-O v\
	-o ${PREFIX}_WGS.vcf\
	$ATAC/WGS/${SAMPLE}_tumor.bam
DIFF=$(($SECONDS - $START))
echo "Time mpileup of WGS: $DIFF seconds"

START=$SECONDS
bcftools mpileup -f $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	-T /broad/medullo/VCF_proteomics_cohort/${SAMPLE}_pair.pass.vcf\
	-a DP,AD\
	-O v\
	-o ${PREFIX}_ATAC.vcf\
	$ATAC/ATAC_bam_hg19/${ABAM}.chr1-22xy.bam
DIFF=$(($SECONDS - $START))
echo "Time mpileup of ATAC: $DIFF seconds"

# Evaluate indels

# Evaluate SNPs
#vcftools --vcf ${PREFIX}_WGS.vcf\
#	--out ${PREFIX}_WGS\
#	--geno-depth
#
#vcftools --vcf ${PREFIX}_ATAC.vcf\
#	--out ${PREFIX}_ATAC\
#	--geno-depth
# I possibly want --counts instead.
# didn't work.
