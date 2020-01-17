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

# Evaluate the quality of the normal samples
source $HOME/.my.bashrc
source activate cnvkit

#cnvkit.py metrics $ATAC/cnvkit/MB*/reference.cnn -s $ATAC/cnvkit/MB*/*.cns > $ATAC/cnvkit/reference_metrics.tsv
cnvkit.py metrics $ATAC/cnvkit/MB*/MB*_control.targetcoverage.cnn > $ATAC/cnvkit/reference_metrics.tsv
cnvkit.py metrics $ATAC/cnvkit/pooled_reference.cnn >> $ATAC/cnvkit/reference_metrics.tsv
conda deactivate

