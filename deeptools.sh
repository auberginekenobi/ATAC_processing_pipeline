#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=4G
# Parallel cores
#$ -pe smp 4
#$ -R y
#$ -binding linear:4


# Runtime
#$ -l h_rt=01:25:00

# Output
#$ -j y

######################
### Dotkit section ###
######################

# This is required to use dotkits inside scripts
#source /broad/software/scripts/useuse

# Use your dotkit
#reuse Bcftools

##################
### Run script ###
##################

# Runs deeptools to generate heatmaps of coverage at specific sites.

# Example usage:
# deeptools.sh -a $ATAC/macs2_hg19/macs2_subcommands/MB018_peaks_wgs_control.bed -b $ATAC/ChIP/enhancers.sorted.GRCh37.bed -o MB018 -A MB018_filterdup.pileup.bw -B H3K27ac/MB18_H3K27Ac_treat_afterfitting_all.bw

source $HOME/.my.bashrc

FLAGS=0
ATAC_BED=
AC_BED=
PREFIX=
ATAC_BW=
AC_BW=
OUTPUT_DIR=$ATAC/deeptools
ATAC_BW_DIR=$ATAC/macs2_hg19/macs2_subcommands/bigwig
AC_BW_DIR=$ATAC/ChIP

# Parse arguments
while [ "$1" != "" ]; do
	case $1 in
		-a | --atac_bed )	shift
					ATAC_BED=""
					while [ $1 != "" ] && [ ${1::1} != "-" ]; do
						ATAC_BED="$ATAC_BED $1"
						shift
					done
					;;
		-b | --ac_bed )		shift
					AC_BED=$1
					shift
					;;
		-o | --output_prefix )	shift
					PREFIX=$1
					shift
					;;
		-A | --atac_bw )	shift
					ATAC_BW=""
					while [ $1 != "" ] && [ ${1::1} != "-" ]; do
						ATAC_BW="$ATAC_BW $ATAC_BW_DIR/$1"
						shift
					done
					;;
		-B | --ac_bw )		shift
					AC_BW=""
					while [ $1 != "" ] && [ ${1::1} != "-" ]; do
						AC_BW="$AC_BW $AC_BW_DIR/$1"
						shift
					done
					;;
		-O | --output_dir )	shift 
					OUTPUT_DIR=$1
					shift
					;;
		-D | --atac_bw_dir )	shift
					ATAC_BW_DIR=$1
					shift
					;;
		-E | --ac_bw_dir )	shift
					AC_BW_DIR=$1
					shift
					;;
	esac
done
echo "ATAC bed files: ${ATAC_BED}"
echo "h3k27ac bed files: ${AC_BED}"
echo "prefix: ${PREFIX}" 
echo "ATAC bw files: $ATAC_BW"
echo "h3k27ac bw files: $AC_BW"

START=$SECONDS

source activate deeptools


OUT=$OUTPUT_DIR/$PREFIX

if [ "$ATAC_BED" != "" ]; then
# Plot matrix by ATAC
computeMatrix reference-point\
	--referencePoint center\
	-a 2000 -b 2000\
	-R $ATAC_BED\
	-S $ATAC_BW $AC_BW\
	-p max\
	-bs 50\
	--missingDataAsZero\
	-o ${OUT}_atac0_mat.gz &&

plotHeatmap -m ${OUT}_atac0_mat.gz\
	-o ${OUT}_atacpeaks_atacsorted.svg\
	--sortUsingSamples 1\
	--heatmapHeight 100\
	--boxAroundHeatmaps no\
	--outFileSortedRegions ${OUT}_atacpeaks_atacsorted.bed\
	--sortUsing max &&

plotHeatmap -m ${OUT}_atac0_mat.gz\
	-o ${OUT}_atacpeaks_acsorted.svg\
	--sortUsingSamples 2\
	--heatmapHeight 100\
	--boxAroundHeatmaps no\
	--outFileSortedRegions ${OUT}_atacpeaks_acsorted.bed\
	--sortUsing max
fi

if [ "$AC_BED" != "" ]; then
# Plot matrix by AC
computeMatrix reference-point\
	--referencePoint center\
	-a 4000 -b 4000\
	-R $AC_BED\
	-S $ATAC_BW $AC_BW\
	-p max\
	-bs 50\
	-o ${OUT}_h3k27ac_mat.gz &&

plotHeatmap -m ${OUT}_h3k27ac_mat.gz\
	-o ${OUT}_h3k27ac_heatmap.svg\
	--sortUsingSamples 2\
	--heatmapHeight 100\
	--boxAroundHeatmaps no\
	--outFileSortedRegions ${OUT}_h3k27ac_sorted.bed\
	--sortUsing max
fi

conda deactivate

DIFF=$(($SECONDS - $START))
echo "deeptools heatmap took $DIFF seconds."
