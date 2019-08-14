#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=6G
# Parallel cores
#$ -pe smp 8
#$ -R y
#$ -binding linear:8

# Runtime
#$ -l h_rt=05:00:00

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

# $1 sample name, eg MB018

bcftools mpileup\
	-f $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	-T /broad/medullo/VCF_proteomics_cohort/${1}_pair.pass.vcf\
	$ATAC/WGS/${1}_tumor.bam
