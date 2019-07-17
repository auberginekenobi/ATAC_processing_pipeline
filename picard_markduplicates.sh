#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=16G

# Runtime
#$ -l h_rt=04:00:00

##################
### Run script ###
##################
# $1 file
# $2 output prefix

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
samtools view -h ${QNAME_SORT_BAM_FILE} | $(which assign_multimappers.py) -k 2 | samtools view -F 1804 -Su /dev/stdin | sambamba sort /dev/stdin -o ${FILT_BAM_FILE}
rm -f ${QNAME_SORT_BAM_FILE}

# ========================
# Mark duplicates
# ======================
FILT_BAM_FILE="${DIR}/${2}.filt.bam"
TMP_FILT_BAM_FILE="${DIR}/${2}.bam"
DUP_FILE_QC="${DIR}/${2}.dup.qc"# QC file
java -Xmx4G -jar $PICARD MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# ============================
# Remove duplicates
# Index final position sorted BAM
#=============================
FINAL_BAM_PREFIX="${DIR}/${2}.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

samtools view -F 1804 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE}

# Index Final BAM file
sambamba index ${FINAL_BAM_FILE}

samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS}

# =============================
# Compute library complexity
# =============================
# sort by position and strand
# Obtain unique count statistics

PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
bedtools bamtobed -i ${FILT_BAM_FILE} | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC}

rm ${FILT_BAM_FILE}

