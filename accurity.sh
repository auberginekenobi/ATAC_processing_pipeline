#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Singularity requires Redhat 7
#$ -l os=RedHat7

# Memory per core
#$ -l h_vmem=6G
# Parallel cores
#$ -pe smp 8
#$ -R y
#$ -binding linear:8

# Runtime
#$ -l h_rt=07:00:00
# MB018 took 5+ hours

# Output
#$ -j y

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
#source /broad/software/scripts/useuse

# Use your dotkit
#reuse Python-3.4

##################
### Run script ###
##################

# $1 Sample id, eg MB018

ATAC=/broad/medullo/ATACseq_fraenkel
OUT=$ATAC/accurity/output/${1}

#mkdir $OUT

singularity exec \
	-B $ATAC/WGS \
	-B /broad/medullo/WGS_fastq_proteomics_cohort/GP_aligned \
	-B /seq/references/Homo_sapiens_assembly19/v1 \
	-B $ATAC/accurity \
	$ATAC/accurity/accurity.simg \
/usr/local/Accurity/main.py \
	-c $ATAC/accurity/configure \
	--nCores 8 \
	-t $ATAC/WGS/${1}_tumor.bam \
	-n $ATAC/WGS/${1}_control.bam \
	-o $OUT \
	--snp_output_dir $OUT \
	-d 1 --clean 1

