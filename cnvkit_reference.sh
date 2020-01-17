#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=12G

# Runtime
#$ -l h_rt=03:00:00
# Parallel cores
# #$ -pe smp 8
# #$ -R y
# #$ -binding linear:8

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

# Generate a pooled reference from all of the .targetcoverage files.

START=$SECONDS

source $HOME/.my.bashrc
source activate cnvkit

ATAC=/broad/medullo/ATACseq_fraenkel

START=$SECONDS
cnvkit.py reference $ATAC/cnvkit/MB*/MB*_control.targetcoverage.cnn\
	--fasta $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	-o $ATAC/cnvkit/pooled_reference.cnn

conda deactivate

END=$SECONDS
DIFF=$(($END - $START))
echo "Reference creation time: $DIFF seconds."

