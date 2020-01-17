# Memory
#$ -l h_vmem=4G

# Runtime
#$ -l h_rt=06:00:00

#$ -j y
#$ -cwd

######################

source /broad/software/scripts/useuse

reuse SRAtoolkit

######################
# Begin work section #
######################
# $1 SRR ID
# $2 output directory
SRR_ID=$1

fastq-dump --outdir $2 --gzip --skip-technical -F --read-filter pass --dumpbase --split-3 --clip -v $SRR_ID
