#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=01:00:00

##################
### Run script ###
##################

# $1 file
# $2 output prefix
# $3 -k param for assign_multimappers.py. Probably 2.
DIR=/broad/medullo/ATACseq_fraenkel/ATAC_bam_hg19
shopt -s expand_aliases
source $HOME/.bashrc
# =============================
Remove unmapped, mate unmapped not primary alignment, reads failing platform 
# ==================
QNAME_SORT_BAM_FILE="${DIR}/${2}.qnmsrt.bam"
FILT_BAM_PREFIX="${DIR}/${2}.filt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"

sambamba sort -n "${DIR}/${1}" -o ${QNAME_SORT_BAM_FILE}
samtools view -h ${QNAME_SORT_BAM_FILE} | $(which assign_multimappers.py) -k $3 | samtools view -F 1804 -Su /dev/stdin | sambamba sort /dev/stdin -o ${FILT_BAM_FILE}
rm -f ${QNAME_SORT_BAM_FILE}
