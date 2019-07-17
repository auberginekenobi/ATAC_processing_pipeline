#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=02:25:00

##################
### Run script ###
##################

# $1 path to input file
# $2 path to output file

source ~/.bashrc
~/bin/hmmcopy_utils/bin/readCounter -c 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y $1 > $2

