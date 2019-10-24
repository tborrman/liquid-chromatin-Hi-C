#!/usr/bin/env python
import argparse

import cooler
import bioframe
import pandas as pd
from copy import copy

import numpy as np
import sys

parser=argparse.ArgumentParser(description='Get cis percent for window range of bin of Hi-C cooler file' +
' for each row of Hi-C matrix')
parser.add_argument('-i', help='input cooler Hi-C file or cooler URI', type=str, required=True)
parser.add_argument('-r', help='window range', type=int, default=6000000)
parser.add_argument('-n', help='maximum percent of NAs allowed in window', type=float, default=10.0)
args=parser.parse_args()


# source of inspiration and copy/paste: https://github.com/mirnylab/cooltools/blob/master/cooltools/coverage.py

def _zero_diags(chunk, n_diags):
    """
    fill n_diags with zeros for a given chunk of pixels
    """
    if n_diags > 0:
        mask = np.abs(chunk['pixels']['bin1_id'] - chunk['pixels']['bin2_id']) < n_diags
        chunk['pixels']['count'][mask] = 0
    return chunk


def _zero_range(chunk, range_bins):
    """
    fill pixels that are separated with more than range_bins apart with zeros
    for a given chunk of pixels
    """
    mask = np.abs(chunk['pixels']['bin1_id'] - chunk['pixels']['bin2_id']) > range_bins
    chunk['pixels']['count'][mask] = 0
    return chunk


def _timesouterproduct(chunk, bias, column_name='balanced'):
    """
    multiply raw pixel counts by the balancing bias
    and return a chunk with an additional column
    """
    pixels = chunk["pixels"]
    chunk['pixels'][column_name] = \
        bias[pixels["bin1_id"]] * bias[pixels["bin2_id"]]* pixels["count"]
        # returning modified chunks with an additional column:
    return chunk

def _get_chunk_coverage_cis(chunk, pixel_weight_key='count'):
    '''
    Compute cisRange "coverage" for a given chunk of pixels.

    Parameters
    ----------
    chunk : dict of dict/pd.DataFrame
        A cooler chunk of pixels produced by the cooler split-apply-combine pipeline.
    pixel_weight_key: str
        The key of a pixel chunk to retrieve pixel weights.

    Returns
    -------
    cov : np.array n_bins
        A numpy array with the cisRange coverage.
    '''
    bins = chunk['bins']
    pixels = chunk['pixels']
    n_bins = len(bins['chrom'])
    cov = np.zeros(n_bins)
    pixel_weights = pixels[pixel_weight_key]

    cis_mask = bins['chrom'][pixels['bin1_id']] == bins['chrom'][pixels['bin2_id']]
    cov += np.bincount(pixels['bin1_id'], weights=pixel_weights*cis_mask, minlength=n_bins)
    cov += np.bincount(pixels['bin2_id'], weights=pixel_weights*cis_mask, minlength=n_bins)

    # skip dividing by 2
    # and beware the main diagonal double counting
    return cov

def main():

    # percent of NAs allowed
    p = args.n/100

    clr = cooler.Cooler(args.i)


    binsize = clr.info['bin-size']
    n_bins = clr.info['nbins']
    grange = args.r/2
    range_bins=int(grange/binsize)
    # chunksize

    # for the cooler-ones ...
    max_NAs = ((tyler_grange/binsize) * NaNsfraction)

    # cooler-specific ones ...
    ignore_diags = 0
    clr_bins = clr.bins()[:]
    bias = clr_bins['weight'].values
    # should we turn NaNs to 0-s ?
    _bias = np.nan_to_num(bias)

    # start with actual calculations,
    # replacing "cp = mf.get_cis_percent_range(f, args.r, p)" ...

    # generate a stream of pixels (split into chunks) and apply various
    chunks = cooler.tools.split(clr, chunksize=int(1e7), map=map, use_lock=False)
    chunks = chunks.pipe(_zero_diags, n_diags=ignore_diags) if ignore_diags>0 else chunks

    cov = (
        chunks
            .pipe(_zero_range, range_bins=range_bins)
            .pipe(_timesouterproduct, bias = _bias)
            .pipe(_get_chunk_coverage_cis, pixel_weight_key="balanced")
            .reduce( np.add, np.zeros(n_bins) )
    )

    # we can do the same stuff + ignore range bins =0
    # and get the double counted main diagonal ...
    # this is very inefficient, but whatever:
    chunks = cooler.tools.split(clr, chunksize=int(1e7), map=map, use_lock=False)
    diag = (
        chunks
            .pipe(_zero_range, range_bins=0)
            .pipe(_timesouterproduct, bias = _bias)
            .pipe(_get_chunk_coverage_cis, pixel_weight_key="balanced")
            .reduce( np.add, np.zeros(n_bins) )
    )

    # afterwards covs-diag/2 should be the same as clr.maxtrix()[:].sum(axis=0) ...

    chrom_pairs = zip(clr.chromnames[:-2],clr.chromnames[1:-1])
    offsets = [(clr.offset(chr1)+range_bins, clr.offset(chr2)-range_bins) for chr1,chr2 in chrom_pairs]

    # some final  modifications we need to do
    # to the cooler-pipeline generated result
    # to make it indistinguishable from Tyler's:

    cp = (cov-diag/2)
    cp[cov==0] = np.nan
    cp[tyler_nans > max_NAs] = np.nan

    # to be rewritten ...
    # Write output
    OUT = open(args.i[:-5] + '_range' + str(args.r/1000000) + 'Mb_cispercent.bedGraph', 'w')
    # Only using 22 autosomes and X
    y_chrom_bin_start =  f['chr_bin_range'][:][23][0]
    for i, b in enumerate(f['bin_positions'][:][:y_chrom_bin_start]):
        if b[0] == 22:
            chrom = 'chrX'
        else:
            chrom = 'chr' + str(b[0]+1)
        OUT.write(chrom+'\t'+str(b[1])+'\t'+str(b[2])+
        '\t'+ str(cp[i])+'\n')
    OUT.close()


if __name__ == '__main__':
    main()
