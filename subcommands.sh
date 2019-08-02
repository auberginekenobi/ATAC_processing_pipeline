#!/bin/bash

#############################
### Default UGER Requests ###
#############################

# Memory
#$ -l h_vmem=32G

# Runtime
#$ -l h_rt=08:00:00

#$ -j y

##################
### Run script ###
##################
source $HOME/.bashrc
source activate py2

cd /broad/medullo/ATACseq_fraenkel/macs2_hg19/macs2_subcommands/

# Params 
# 1 Sample name, eg MB248
# 2 [optional] if the ATAC bam file is named something else, add it as a second argument.
# Usage:
# subcommands.sh MB106 MB106c

sample=$1
atac=${2:-$1}

# 1 Filter duplicates
macs2 filterdup -i ../../ATAC_bam_hg19/${atac}.chr1-22xy.bam --keep-dup=all -o ${atac}_temp_filterdup.bed 2> ${atac}_filterdup_stats.txt
E1="$(tr -d '\n' < ${atac}_filterdup_stats.txt)"; wait #strip newlines from log
E1="$(expr match "$E1" '.*\(tags in alignment file: [0-9]*\)')" # match correct line
E1="$(expr match "$E1" '[^0-9]*\([0-9]*\)')" # get value
echo "reads in ATACseq file: $E1"; wait

macs2 filterdup -i ../../WGS/${sample}_control.bam --keep-dup=1 -o ${sample}_temp_control_filterdup.bed 2> ${sample}_temp_control.txt
E2="$(tr -d '\n' < ${sample}_temp_control.txt)"; wait
E2="$(expr match "$E2" '.*\(tags after filtering in alignment file: [0-9]*\)')"
E2="$(expr match "$E2" '[^0-9]*\([0-9]*\)')"
echo "reads in WGS file: $E2"; wait

# 2 

# 3 Extend ChIP sample to get ChIP coverage track
macs2 pileup -i ${atac}_temp_filterdup.bed -o ${atac}_filterdup.pileup.bdg -B --extsize 100 -f BED
echo "Generated coverage track"

# 4 Build local bias tracks

# 4.1 Size d background (200)
macs2 pileup -i ${sample}_temp_control_filterdup.bed -B --extsize 100 -o ${sample}_temp_d_bg.bdg -f BED
echo "Generated d background"

# 4.2 

# 4.3 The llocal background
macs2 pileup -i ${atac}_temp_filterdup.bed -B --extsize 5000 -o ${atac}_temp_10k_bg.bdg -f BED
macs2 bdgopt -i ${atac}_temp_10k_bg.bdg -m multiply -p 0.02 -o ${atac}_temp_10k_bg_norm.bdg
rm ${atac}_temp_10k_bg.bdg
echo "Generated 10k background"

# 4.4 The whole genome background
# The WG background can be calcualted as n_control_reads*frag_length/genome_size
WGB="$(echo "scale=5;${E2}*200/2700000000" | bc)"; wait
echo "Generated WGB background: $WGB"

# 4.5 Combine and generate maximum background noise
macs2 bdgopt -i ${sample}_temp_d_bg.bdg -m max -p $WGB -o ${sample}_temp_local_bias_1.bdg
echo "Combining backgrounds..."

# 5. Scale ChIP and Control to same sequencing depth
# -p = ATAC reads / control reads
SCALE_FACTOR="$(echo "scale=5;${E1}/${E2}" | bc)"; wait
echo "Scale factor: $SCALE_FACTOR"
macs2 bdgopt -i ${sample}_temp_local_bias_1.bdg -m multiply -p ${SCALE_FACTOR} -o ${sample}_temp_local_bias_2.bdg
echo "Scaling backgrounds..."

# 4.5 Combine final background noises
macs2 bdgcmp -m max -t ${atac}_temp_10k_bg_norm.bdg -c ${sample}_temp_local_bias_2.bdg -o ${atac}_local_lambda_atac_llocal.bdg
echo "Done combining backgrounds"

#6. Compare ChIP and local lambda to get the scores in pvalue or qvalue 
macs2 bdgcmp -t ${atac}_filterdup.pileup.bdg -c ${atac}_local_lambda_atac_llocal.bdg -m qpois -o ${atac}_temp_qvalue.bdg
echo "Generated q-value track"

#7. Call peaks on score track
macs2 bdgpeakcall -i ${atac}_temp_qvalue.bdg -c 2 -l 200 -g 50 -o ${atac}_peaks_wgs_control.bed
echo "Cleaning up..."

rm ${atac}_temp_*
rm ${sample}_temp_*
echo "Done!"
conda deactivate
