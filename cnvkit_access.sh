#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=00:30:00

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

# Call before running cnvkit.sh to generate the access file.

source $HOME/.my.bashrc
source activate cnvkit

ATAC=/broad/medullo/ATACseq_fraenkel
PREFIX=$ATAC/cnvkit/Homo_sapiens_assembly19

cnvkit.py access $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	-o $PREFIX.temp.bed

# Keep only chromosomes 1-22 XY
sed -n '/^[0-9,X,Y]/Ip' $PREFIX.temp.bed > $PREFIX.access.bed
rm $PREFIX.temp.bed
conda deactivate

