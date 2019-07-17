#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=04:00:00

##################
### Run script ###
##################
# $1 trimmomatic fastq file: MB274_SHH_SHHa_BATCH3-170809Fra_D17-6964.fastq
# $2 output sam file: MB***.sam

DIR=/broad/medullo/ATACseq_fraenkel

source $HOME/.bashrc

cd $DIR
echo $2
bowtie2 -x genomes/hg19 --phred33 -U trimmomatic/$1 -S ATAC_bam_hg19/$2
