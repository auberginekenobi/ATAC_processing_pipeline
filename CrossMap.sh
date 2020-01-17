#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=00:45:00

# Output
#$ -j y

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
#source /broad/software/scripts/useuse

# Use your dotkit
#reuse Bcftools

##################
### Run script ###
###################

# UGER wrapper script around crossmap
# Usage: qsub CrossMap.sh bigwig $ATAC/genomes/hg19ToGRCh37.over.chain.gz $ATAC/ChIP/H3K27ac/MB95_H3K27Ac_treat_afterfiting_all.bw $ATAC/ChIP/H3K27ac/MB95_H3K27Ac_GRCh37

source .my.bashrc
source activate CrossMap

START=$SECONDS
CrossMap.py $@
DIFF=$(($SECONDS - $START))
echo "CrossMap took $DIFF seconds."
conda deactivate
