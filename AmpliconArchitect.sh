#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Singularity requires Redhat 7
#$ -l os=RedHat7

# Memory per core
#$ -l h_vmem=12G
# Parallel cores

# Runtime
#$ -l h_rt=35:00:00
# MB018 took 5+ hours

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

# Run after merge_alt.py; assumes MBxxx_merge.dat
# $1 Sample id, eg MB018

source $HOME/.my.bashrc
ATAC=/broad/medullo/ATACseq_fraenkel
AA=/broad/medullo/ATACseq_fraenkel/AmpliconArchitect/AmpliconArchitect/src/AmpliconArchitect.py

input_bam=$ATAC/WGS/${1}_tumor.bam
bed_file=$ATAC/AmpliconArchitect/${1}/${1}_merge.dat
output_prefix=$ATAC/AmpliconArchitect/${1}/${1}_AA
rm ${output_prefix}*

singularity exec \
	-B $ATAC/WGS \
	-B $ATAC/AmpliconArchitect \
	-B $HOME/mosek \
	-B /broad/medullo/WGS_fastq_proteomics_cohort/GP_aligned \
	$ATAC/AmpliconArchitect/AmpliconArchitect.simg \
$AA --bam ${input_bam} --bed ${bed_file} --out ${output_prefix} \
	--ref GRCh37 
#	--downsample -1
