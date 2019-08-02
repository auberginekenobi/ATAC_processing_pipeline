#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=00:45:00

##################
### Run script ###
##################

# $1 file prefix

DIR=/broad/medullo/ATACseq_fraenkel/ATAC_bam_hg19
source $HOME/.bashrc

cd $DIR
echo $1

#samtools index "${1}.nodup.bam"

## Check if index is present. If not, create it:
if [[ ! -e ${1}.nodup.bam.bai ]];
  then
  echo '[INFO]: File does not seem to be indexed. Indexing now:'
  samtools index ${1}.nodup.bam
  fi

## Calculate %mtDNA:
mtReads=$(samtools idxstats "${1}.nodup.bam" | grep 'chrM\|MT' | cut -f 3)
totalReads=$(samtools idxstats "${1}.nodup.bam" | awk '{SUM += $3} END {print SUM}')

echo '==> mtDNA Content:' $(bc <<< "scale=2;100*$mtReads/$totalReads")'%'

## Remove contigs, mitochondria
samtools view -bh "${1}.nodup.bam" $(cat $HOME/scripts/hg19_chr.tsv) > "${1}.chr1-22xy.bam" 
samtools index "${1}.chr1-22xy.bam"
samtools flagstat "${1}.chr1-22xy.bam" > "${1}.chr1-22xy.bam.flagstat.qc"
echo "Done."
 
