#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=6G

# Runtime
#$ -l h_rt=02:00:00

#$ -j y
#$ -cwd
##################
### Run script ###
##################
source $HOME/.bashrc
source activate py2

# $1 full path to input bam
# $2 output prefix (optional)

genomesize=2700000000
DIR="$( dirname "$1" )"
base="$( basename "$1" .chr1-22xy.bam )"
NAME=${2:-$base}

macs2 callpeak \
	-t $1 \
	--name $NAME \
	--format BAM \
	--gsize $genomesize \
	--slocal 1000 \
	--llocal 10000 \
	--call-summits \
	--keep-dup 'all' \
	--bdg \
	--qvalue '0.01' \
	--nomodel \
	--extsize '200' \
	--shift '-100' 
echo "Done!"
conda deactivate
