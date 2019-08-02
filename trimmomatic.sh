#!/bin/bash
# Memory
#$ -l h_vmem=2G

# Runtime
#$ -l h_rt=01:00:00

#$ -j y

######################
# Begin work section #
######################
# $1 full path to fastq file input
# $2 output file prefix, eg MB***
# output goes to DIR/trimmomatic/$2.fastq

DIR=/broad/medullo/ATACseq_fraenkel

source $HOME/.bashrc

java -jar $TRIMMOMATIC SE $1 $DIR/trimmomatic/$2.fastq LEADING:15 TRAILING:15
