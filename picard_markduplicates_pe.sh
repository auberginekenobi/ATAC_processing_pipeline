#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=16G

# Runtime
#$ -l h_rt=04:00:00

#$ -j y

##################

source /broad/software/scripts/useuse

reuse BEDTools

##################
### Run script ###
##################
# $1 file

DIR="$( dirname "$1" )"
base="$( basename "$1" .sorted.bam )"
shopt -s expand_aliases
source $HOME/.bashrc
source activate pybedtools

echo $DIR
echo $base

function fail {
	echo "$@" >&2
	exit 1
}
set -o pipefail

# =============================
# Remove unmapped, mate unmapped not primary alignment, reads failing platform
# =============================
FILT_BAM_PREFIX="${DIR}/${base}.filt"
FILT_BAM_FILE="${FILT_BAM_PREFIX}.bam"
TMP_FILT_BAM_PREFIX="${FILT_BAM_PREFIX}.tmp"
TMP_FILT_BAM_FILE="${TMP_FILT_BAM_PREFIX}.bam"

samtools view -h -F 524 -f 2 -u $1 | samtools sort -n /dev/stdin -o ${TMP_FILT_BAM_FILE} || 
	fail "Is samtools working?"

TMP_FILT_FIXMATE_BAM_FILE="${TMP_FILT_BAM_PREFIX}.fixmate.bam"
samtools view -h ${TMP_FILT_BAM_FILE} | 
	$(which assign_multimappers.py) -k 2 --paired-end | 
	samtools fixmate -r /dev/stdin ${TMP_FILT_FIXMATE_BAM_FILE} ||
	fail "Samtools fixmate failed"

# =============================
# Remove orphan reads (pair was removed)
# and read pairs mapping to different chromosomes
# Obtain position sorted BAM
# ============================= 
samtools view -F 1804 -f 2 -u ${TMP_FILT_FIXMATE_BAM_FILE} | samtools sort /dev/stdin -o ${FILT_BAM_FILE} ||
	fail "Sort and filter failed"
rm ${TMP_FILT_FIXMATE_BAM_FILE}
rm ${TMP_FILT_BAM_FILE} 

# ========================
# Mark duplicates
# ======================
DUP_FILE_QC="${FILT_BAM_PREFIX}.dup.qc" # QC file
#samtools markdup -s ${FILT_BAM_FILE} ${TMP_FILT_BAM_FILE} > ${DUP_FILE_QC} 2>&1

java -Xmx4G -jar $PICARD MarkDuplicates INPUT=${FILT_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false ||
	fail "picard markduplicates failed."
mv ${TMP_FILT_BAM_FILE} ${FILT_BAM_FILE}

# ============================
# Remove duplicates
# Index final position sorted BAM
#=============================
FINAL_BAM_PREFIX="${DIR}/${base}.nodup"
FINAL_BAM_FILE="${FINAL_BAM_PREFIX}.bam" # To be stored
FINAL_BAM_INDEX_FILE="${FINAL_BAM_FILE}.bai" # To be stored
FINAL_BAM_FILE_MAPSTATS="${FINAL_BAM_PREFIX}.flagstat.qc" # QC file

# Final bam contains only reads with flags 1804 unset.
samtools view -F 1804 -f 2 -b ${FILT_BAM_FILE} > ${FINAL_BAM_FILE} || fail "Samtools 1804 failed"

# Index Final BAM file
samtools index ${FINAL_BAM_FILE} || fail "Samtools indexing failed"

samtools flagstat ${FINAL_BAM_FILE} > ${FINAL_BAM_FILE_MAPSTATS} || fail "Samtools flagstat failed"

# =============================
# Compute library complexity
# =============================
# sort by position and strand
# Obtain unique count statistics

PBC_FILE_QC="${FINAL_BAM_PREFIX}.pbc.qc"

# PBC File output
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] PBC2=OnePair/TwoPair
samtools sort -n ${FILT_BAM_FILE} -o ${TMP_FILT_BAM_FILE}
bedtools bamtobed -bedpe -i ${TMP_FILT_BAM_FILE} | 
	awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' | 
	grep -v 'chrM\|MT' | sort | uniq -c | 
	awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > ${PBC_FILE_QC} ||
	fail "QC metrics failed."

rm ${FILT_BAM_FILE}
rm ${FILT_BAM_FILE}.bai
rm ${TMP_FILT_BAM_FILE}
conda deactivate
echo "Done!"
