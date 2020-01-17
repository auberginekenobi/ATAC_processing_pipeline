#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=6G
# Parallel cores
#$ -pe smp 4
#$ -R y
#$ -binding linear:4

# Runtime
#$ -l h_rt=00:45:00

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

# $1 Sample ID, eg MB018
# $2 [optional] "tumor" or "control". Default "tumor"
SAMPLE=$1
SUF=${2:-"tumor"}

START=$SECONDS

source $HOME/.my.bashrc
source activate cnvkit

ATAC=/broad/medullo/ATACseq_fraenkel
OUT=$ATAC/cnvkit/${1}

cnvkit.py fix $OUT/${SAMPLE}_${SUF}.targetcoverage.cnn\
	$OUT/${SAMPLE}_${SUF}.antitargetcoverage.cnn\
	$ATAC/cnvkit/pooled_reference.cnn\
	-o $OUT/${SAMPLE}_temp.cnr

# Keep only chromosomes 1-22 XY
head -n 1 $OUT/${SAMPLE}_temp.cnr > $OUT/${SAMPLE}_${SUF}.cnr
sed -n '/^[0-9,X,Y]/Ip' $OUT/${SAMPLE}_temp.cnr >> $OUT/${SAMPLE}_${SUF}.cnr
rm $OUT/${SAMPLE}_temp.cnr

cnvkit.py segment $OUT/${SAMPLE}_${SUF}.cnr -o $OUT/${SAMPLE}_${SUF}.cns\
	--drop-low-coverage\
	--drop-outliers 15\
	-t 1e-7\
	-p 4


cnvkit.py scatter $OUT/${SAMPLE}_${SUF}.cnr -s $OUT/${SAMPLE}_${SUF}.cns -o $OUT/${SAMPLE}_${SUF}-scatter.png
#cnvkit.py scatter $OUT/${1}_tumor.cnr -s $OUT/${1}_tumor.cns -o $OUT/${1}_tumor-scatter.pdf
# diagram can only generate pdf
# This one takes a long time with the .cnr
cnvkit.py diagram -s $OUT/${SAMPLE}_${SUF}.cns -o $OUT/${SAMPLE}_${SUF}-diagram.pdf

conda deactivate

END=$SECONDS
DIFF=$(($END - $START))
echo "Processing time: $DIFF seconds."
