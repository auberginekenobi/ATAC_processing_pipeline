#!/bin/bash
# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=02:00:00

#$ -j y

######################
# Begin work section #
######################
# $1 SAF file to count
# $2 list of files in the form '"file1 file2 ... filen"'
source .bashrc

DIR=/broad/medullo/ATACseq_fraenkel
BAM=$DIR/ATAC_bam_hg19
INPUT=$DIR/featureCounts
OUTPUT=$DIR/featureCounts
SAF=$INPUT/$1
cd $BAM


echo -n $2 | xargs featureCounts -a $SAF -o $OUTPUT/narrow_peak_counts.raw.txt -F SAF -T 4

