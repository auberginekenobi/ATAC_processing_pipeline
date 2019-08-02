#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=00:30:00

#$ -j y

##################
### Run script ###
##################

# $1 bam file
# $2 output prefix
# $3 read length
# Usage:
# tss_enrichment.sh MB106c.chr1-22xy.bam MB106c  45

source $HOME/.bashrc
source activate pybedtools
$HOME/scripts/tss_enrichment.py -b $1 -o $2 -l $3
conda deactivate
