#!/bin/bash
# Memory
#$ -l h_vmem=2G

# Runtime
#$ -l h_rt=01:00:00

#$ -cwd
#$ -j y

######################
# Begin work section #
######################
# $1 full path to fastq file input
# $2 full path to fastq read file 2
# $3 output file prefix, eg MB***

source $HOME/.bashrc

path=$(dirname "$3")
if [ $path == '.' ]; then
	path=$(dirname "$1")
fi
out=$(basename "$3" .fastq)

java -jar $TRIMMOMATIC PE $1 $2 \
	$path/${out}_paired_1.fastq \
	$path/${out}_unpaired_1.fastq \
	$path/${out}_paired_2.fastq \
	$path/${out}_unpaired_2.fastq \
	LEADING:15 TRAILING:15
