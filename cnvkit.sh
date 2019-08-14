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
#$ -l h_rt=04:00:00

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

# $1 Sample ID, eg MB018

source $HOME/.my.bashrc
source activate cnvkit

ATAC=/broad/medullo/ATACseq_fraenkel

cnvkit.py batch $ATAC/WGS/${1}_tumor.bam\
	-n $ATAC/WGS/${1}_control.bam\
	-m wgs\
	-f $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	-d $ATAC/cnvkit/${1}\
	--scatter --diagram\
	-p 8\
	--target-avg-size 1000\
	-g $ATAC/cnvkit/Homo_sapiens_assembly19.access.bed
#	-t $ATAC/cnvkit/Homo_sapiens_assembly19.target.1000.bed
#	-a $ATAC/cnvkit/Homo_sapiens_assembly19.antitarget.1000.bed

conda deactivate

