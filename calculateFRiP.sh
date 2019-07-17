#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=00:30:00

##################
### Run script ###
##################

# $1 peak file
# $2 bam file

DIR=/broad/medullo/ATACseq_fraenkel
BAM="${DIR}/ATAC_bam_hg19"
PEAKS="${DIR}/macs2_hg19/macs2_subcommands"
FRIP="${DIR}/macs2_hg19/macs2_subcommands/frip"
source $HOME/.bashrc

cd $DIR
echo $1



awk 'BEGIN{OFS="\t"} FNR>1 {print $1":"$2+1"-"$3, $1, $2+1, $3, "."}' "${PEAKS}/${1}" > "${PEAKS}/${1}.saf"
## run featureCounts
featureCounts -p -a "${PEAKS}/${1}.saf" -F SAF -o "${FRIP}/${1}" "${BAM}/${2}" -T 4
#rm "${PEAKS}/${1}.saf"
