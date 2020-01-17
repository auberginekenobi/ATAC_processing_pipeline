#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory per core
#$ -l h_vmem=4G

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
	echo "	bigWigCrossMapToBedGraph.sh - Convert bigwig to bedgraph in the desired assembly."
	echo "USAGE"
	echo "	bigWigCrossMapToBedGraph [options] bigwig"
	exit 1
}

CHAIN=
while getopts ":o:c:h" opt; do
	case $opt in
		o)
		OUTPUT=$OPTARG
		;;
		c)
		CHAIN=$OPTARG
		;;
		\?)
		echo -e "Invalid option: -$OPTARG\n\n" >&2
		usage >&2
		;;
		:)
		echo -e "Option -$OPTARG requires an argument\n\n" >&2
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


path=$(dirname -- "$INPUT")
ext=${INPUT##*.}
base=$(basename -- "$INPUT" ".$ext")
if [ "$OUTPUT" == "" ]; then
	OUTPUT=$path/$base.bdg
fi

if [ "$CHAIN" == "" ]; then
	TEMP=$OUTPUT
else
	[ -e $CHAIN ] || fail "could not find $CHAIN"
	TEMP=$path/$base.temp.bdg
fi

echo "INPUT FILE: $INPUT"
echo "OUTPUT FILE: $OUTPUT"
echo "CHAIN FILE: $CHAIN"

[ -e $INPUT ] || fail "Could not find $INPUT"

# Convert to bedgraph
echo "Converting to bedgraph..."
START=$SECONDS
bigWigToBedGraph $INPUT $TEMP \
	|| fail "bigWigToBedGraph failed"
DIFF=$(($SECONDS - $START))
echo "Done. bigWigToBedGraph took $DIFF seconds."

# Chain to desired assembly
if [ "$CHAIN" != "" ]; then
	echo "Converting between assemblies using $CHAIN ..."
	START=$SECONDS
	source activate CrossMap
	CrossMap.py bed $CHAIN $TEMP $OUTPUT \
		|| fail "CrossMap $CHAIN $TEMP $OUTPUT failed"
	rm $TEMP
	conda deactivate
	DIFF=$(($SECONDS - $START))
	echo "Done. CrossMap took $DIFF seconds."
fi

