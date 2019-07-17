#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=2G

# Runtime
#$ -l h_rt=01:30:00

##################
### Run script ###
##################

# $1 file
DIR=/broad/medullo/ATACseq_fraenkel/ATAC_bam_hg19

source $HOME/.bashrc

cd $DIR
echo $1
filename=$(basename $1 .sam)
samtools view -Shb $filename.sam > $filename.bam
echo "done viewing"
samtools sort $filename.bam > $filename.sorted.bam
echo "done sorting"
rm $filename.bam
