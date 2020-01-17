#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=6G

# Runtime
#$ -l h_rt=05:00:00

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

# Generate a pooled reference from all of the .targetcoverage files.

START=$SECONDS

#source $HOME/.my.bashrc
#source activate cnvkit

ATAC=/broad/medullo/ATACseq_fraenkel

for file in $ATAC/WGS/MB*_control.bam.bai; do
	ST=$SECONDS
	file="$(basename -- $file)"
	SAMPLE=${file:0:5}
	
	# cnvkit.py autobin -m wgs

	#cnvkit.py coverage $ATAC/WGS/${SAMPLE}_control.bam\
	#	$ATAC/cnvkit/Homo_sapiens_assembly19.target.1000.bed\
	#	-p 8\
	#	-o $ATAC/cnvkit/$SAMPLE/${SAMPLE}_control.targetcoverage.cnn
	DIFF=$(($SECONDS - $ST))
	echo "cnvkit coverage of $SAMPLE took $DIFF seconds."
done

ST=$SECONDS
#cnvkit.py reference $ATAC/cnvkit/MB*/MB*_control.targetcoverage.cnn\
	#--fasta $ATAC/genomes/Homo_sapiens_assembly19.fasta\
	#-o $ATAC/cnvkit/pooled_reference.cnn
DIFF=$(($SECONDS - $ST))
echo "cnvkit reference took $DIFF seconds."

#conda deactivate

END=$SECONDS
DIFF=$(($END - $START))
echo "Reference creation time: $DIFF seconds."

