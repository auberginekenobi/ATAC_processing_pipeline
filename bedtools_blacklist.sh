#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=5G

# Runtime
#$ -l h_rt=00:25:00

##################
### Run script ###
##################

# $1 narrowPeak file

DIR=/broad/medullo/ATACseq_fraenkel/macs2_hg19/q0.01
SOURCE=$DIR/narrowPeak
DEST=$DIR/whitelist_bed
source $HOME/.bashrc
BLACKLIST=$HOME/scripts/wgEncodeDacMapabilityConsensusExcludable.bed

bedtools intersect -v -a $SOURCE/$1 -b $BLACKLIST > $DEST/$1
