#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=08:00:00

##################
### Run script ###
##################

# $1 treatment file 
# $2 control file

DIR=/broad/medullo/ATACseq_fraenkel
SOURCE="${DIR}/ATAC_bam_hg19"
DEST="${DIR}/macs2_hg19"
WGS="${DIR}/../WGS_fastq_proteomics_cohort/GP_aligned/${2}_tumor/${2}_tumor.bam"
source $HOME/.bashrc
use MACS2

cd $DIR
echo $1

macs2 callpeak -t "${SOURCE}/${1}" --control $WGS --name "${2}" --format BAM --gsize '2700000000' --call-summits --slocal 1000 --llocal 10000 --keep-dup '1' --bdg --qvalue '0.01' --nomodel --extsize '200' --shift '-100' --outdir 'macs2_hg19'
conda deactivate
