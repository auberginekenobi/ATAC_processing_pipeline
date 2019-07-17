#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=00:30:00

##################
### Run script ###
##################

# $1 file prefix

DIR=/broad/medullo/ATACseq_fraenkel/ATAC_bam_hg19
source $HOME/.bashrc

cd $DIR
echo $1

samtools index "${1}.nodup.bam"
samtools view -bh "${1}.nodup.bam" $(cat $HOME/scripts/hg19_chr.tsv) > "${1}.chr1-22xy.bam" 
echo "Done."
 
