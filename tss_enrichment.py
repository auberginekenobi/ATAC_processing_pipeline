#!/home/unix/ochapman/miniconda3/envs/pybedtools/bin/python
import matplotlib
matplotlib.use('Agg')

import pybedtools
import metaseq
import logging
import argparse
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import mlab

def main():
    parser = argparse.ArgumentParser(description="This script calculates TSS enrichment given an ATAC-seq bam file.")
    parser.add_argument('-b','--bam', help='Input bam file',required=True)
    parser.add_argument('-o','--output',help='Prefix for output files',required=True)
    parser.add_argument('-g','--genome',help="Genome assembly. Allowable values: 'GRCh37','b37','hg19'. Default 'GRCh37'",required=False,default='GRCh37')
    parser.add_argument('-l','--read_length',help='ATAC-seq read lengths',required=True)
    # TODO add optional arguments
    args = parser.parse_args()
    
    # Assume the Broad bam directory if full path not given
    if not '/' in args.bam:
        args.bam = '/broad/medullo/ATACseq_fraenkel/ATAC_bam_hg19/'+args.bam
    if not '/' in args.output:
        args.output = '/broad/medullo/ATACseq_fraenkel/ATAC_bam_hg19/'+args.output
    
    # Parse genomes
    if args.genome == 'GRCh37' or args.genome == 'b37':
        tss = '/broad/medullo/ATACseq_fraenkel/genomes/gencode.v19.annotation.b37.TSS.bed'
        chromsizes = '/broad/medullo/ATACseq_fraenkel/genomes/human.b37.genome'
    elif args.genome == 'hg19':
        tss = '/broad/medullo/ATACseq_fraenkel/genomes/gencode.v19.annotation.TSS.bed'
        chromsizes = '/broad/medullo/ATACseq_fraenkel/genomes/human.hg19.genome'
    else:
        raise ValueError('Illegal genome argument.')

    args.read_length = int(args.read_length)
    
    print('''Running TSS enrichment with the following arguments:
    bam={}
    output={}
    tss={}
    chromsizes={}
    read_length={}
    '''.format(args.bam,args.output,tss,chromsizes,args.read_length))
    (_,_,tss_point_val)=make_tss_plot(args.bam, tss, args.output, chromsizes, args.read_length)
    print('TSS score({}): {}'.format(args.output,tss_point_val))

def make_tss_plot(bam_file, tss, prefix, chromsizes, read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    '''
    Take bootstraps, generate tss plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating tss plot...')
    tss_plot_file = '{0}_tss-enrich.pdf'.format(prefix)
    tss_plot_large_file = '{0}_large_tss-enrich.pdf'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load the bam file
    bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                          processes=processes, stranded=True)

    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    # Note the middle high point (TSS)
    tss_point_val = max(bam_array.mean(axis=0))

    ax.set_xlabel('Distance from TSS (bp)')
    ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(tss_plot_file)

    # Print a more complicated plot with lots of info

    # Find a safe upper percentile - we can't use X if the Xth percentile is 0
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    plt.rcParams['font.size'] = 8
    fig = metaseq.plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(tss_plot_large_file)

    return tss_plot_file, tss_plot_large_file, tss_point_val

main()
