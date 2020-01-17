#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=6G

# Runtime
#$ -l h_rt=01:15:00
# Parallel cores
#$ -pe smp 8
#$ -R y
#$ -binding linear:8

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

# Generate a coverage file from access and targets (access and autobin)

# $1 bam file in $ATAC/WGS; should have an index. Eg. MB018_tumor.bam
FILE=$1
ROOT="$(basename -s .bam $FILE)"
SAMPLE=${ROOT:0:5}
START=$SECONDS

source $HOME/.my.bashrc
source activate cnvkit

ATAC=/broad/medullo/ATACseq_fraenkel

mkdir -p $ATAC/cnvkit/$SAMPLE
cnvkit.py coverage $ATAC/WGS/${ROOT}.bam\
	$ATAC/cnvkit/Homo_sapiens_assembly19.target.1000.bed\
	-p 8\
	-o $ATAC/cnvkit/$SAMPLE/$ROOT.targetcoverage.cnn

# No antitargets in WGS, but cnvkit seems to require it
cp $ATAC/cnvkit/antitargetcoverage.cnn $ATAC/cnvkit/$SAMPLE/$ROOT.antitargetcoverage.cnn

DIFF=$(($SECONDS - $START))
echo "cnvkit coverage of $SAMPLE took $DIFF seconds."

conda deactivate

