#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=4G
# Parallel cores
#$ -pe smp 4
#$ -R y
#$ -binding linear:4

# Runtime
#$ -l h_rt=04:00:00

#$ -cwd
#$ -j y

##################
### Run script ###
##################
# $1 trimmomatic fastq prefix: everything before _{12}.fastq

base=$( basename $1 )
DIR=$( dirname $1 )

source $HOME/.bashrc

cd $DIR
echo $base
bowtie2 -X 2000 -x /broad/medullo/ATACseq_fraenkel/genomes/hg19 --phred33 \
	-1 ${base}_1.fastq -2 ${base}_2.fastq \
	-S ${base}.sam --threads 4
echo "Done!"
