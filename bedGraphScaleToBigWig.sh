OB#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=8G

# Runtime
#$ -l h_rt=00:15:00

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
###################

# Convert an input file to bedgraph and CrossMap to assembly of choice.

source $HOME/.my.bashrc
set -o pipefail

function fail {
	echo "$@" >&2
	exit 1
}

function usage {
	echo "NAME"
	echo "	bedGraphScaleToBigWig.sh - Scale a bedgraph by the desired amount and return a bigwig.."
	echo "USAGE"
	echo "	bedGraphScaleToBigWig.sh [options] bedgraph genome"
	echo "OPTIONS"
	echo "	-s: Scale factor. Default 1."
	echo "	-o: Output file. Default Input.scale.bw"
	echo "	-r: number of mapped reads. Will calculate scale factor from this assuming RPM. Do not use with -s."
	exit 1
}

SF=1
OUTPUT=
while getopts ":r:o:s:h" opt; do
	case $opt in
		o)
		OUTPUT=$OPTARG
		;;
		s)
		SF=$OPTARG
		;;
		r)
		SF=$(echo "scale=5;20000000 / $OPTARG" | bc -l)
		;;
		\?)
		echo -e "Invalid option: $OPTARG\n\n" >&2
		usage >&2
		;;
		:)
		echo -e "Option -$opt requires an argument\n\n" >&2
		usage >&2
		;;
		h)
		usage
		exit 0
		;;
	esac
done

shift $((OPTIND - 1))
INPUT=${1:-none}
[ "$INPUT" != "none" ] || usage >&2

shift
GENOME=${1:-none}
[ "$GENOME" != "none" ] || usage >&2

path=$(dirname -- "$INPUT")
ext=${INPUT##*.}
base=$(basename -- "$INPUT" ".$ext")
if [ "$OUTPUT" == "" ]; then
	OUTPUT=$path/$base.scale.bw
fi

echo "INPUT FILE: $INPUT"
echo "OUTPUT FILE: $OUTPUT"
echo "SCALE FACTOR: $SF"

[ -e $INPUT ] || fail "Could not find $INPUT"
[ -e $GENOME ] || fail "Could not find $GENOME"

echo "Scaling and sorting..."
START=$SECONDS
TEMP=$path/$base.scale.bdg
export LC_COLLATE=C
if [ $SF != 1 ]; then
	awk -v f=$SF 'BEGIN{OFS="\t"}{$4=f*$4; print}' $INPUT \
		| sort -S 7G -k1,1 -k2,2n > $TEMP \
		|| fail "Scaling or sorting failed"
#	source activate py2 || fail "MACS2 not found"
#	macs2 bdgopt -i $INPUT -m multiply -p $SF -o $TEMP || fail "bdgopt failed"
#	conda deactivate
else
	sort -S 7G -k1,1 -k2,2n $INPUT > $TEMP \
		|| fail "Sorting failed"
fi
DIFF=$(($SECONDS - $START))
echo "Done. Scaling and sorting took $DIFF seconds."


echo "Converting to bigwig..."
START=$SECONDS
bedGraphToBigWig $TEMP $GENOME $OUTPUT || fail "bedGraphToBigWig $TEMP $GENOME $OUTPUT failed"
rm $TEMP
DIFF=$(($SECONDS - $START))
echo "Done. bedGraphToBigWig took $DIFF seconds."
