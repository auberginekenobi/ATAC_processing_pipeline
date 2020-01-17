#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=12G

# Runtime
#$ -l h_rt=01:15:00

# Output
#$ -j y

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
#source /broad/software/scripts/useuse

# Use your dotkit
#reuse BEDTools

##################
### Run script ###
##################

source $HOME/.my.bashrc
source activate py2

# Params
# 1 Sample name, eg MB248

SAMPLE=$1
DIR=/broad/medullo/ATACseq_fraenkel/macs2_hg19/macs2_subcommands

START=$SECONDS
macs2 bdgcmp -t $DIR/${SAMPLE}_filterdup.pileup.bdg -c $DIR/${SAMPLE}_local_lambda_atac_llocal.bdg\
	-o $DIR/${SAMPLE}_FE.bdg\
	-m FE
DIFF=$(($SECONDS - $START))
echo "Time to calculate fold enrichment track: $DIFF seconds."

conda deactivate
