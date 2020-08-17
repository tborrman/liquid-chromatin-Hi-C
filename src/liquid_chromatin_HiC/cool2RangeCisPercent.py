#!/usr/bin/env python
import argparse

import cooler
import bioframe
import pandas as pd
from copy import copy

import numpy as np
import sys


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

def _get_mask(n_bins,
            max_nan,
            bias,
            chromnames,
            chrom_offsets,
            cis_range_size):
    """
    put calculations of the mask into a function

    n_bins - number of bins - size of the mask ultimately
    max_nan - number of NaNs allowed within cis_range on each row
    bias - balancing weights with NaNs for filtered out rows/cols
    chromnames - list of chromosome names
    chrom_offsets - a function that returns chromosome offests
    cis_range_size - a half-width cisRange size

    Return:
    --------------------------------
    mask - boolear array with n_bins valuess
    """

    # we will create a mask to fill cisRange
    # with NaNs where mask is True:

    # initialize with everyhting maked:
    mask = np.ones(n_bins,dtype=np.bool)

    # (1) mask grange of each chrom's start and end
    chrom_pairs = zip(chromnames[:-2],chromnames[1:-1])
    chrom_offsets = [ (chrom_offsets(chr1)+cis_range_size, \
                       chrom_offsets(chr2)-cis_range_size) \
                    for chr1,chr2 in chrom_pairs]
    # unmask the middle of each chromosome:
    for _1,_2 in chrom_offsets:
        mask[_1:_2] = False
    
    # (2) rows/cols filtered by balancing should be masked
    mask[np.isnan(bias)] = True

    # (3) mask rows/cols with too many NaNs in cisRange:
    # count # of NaNs in a given cisRange for a given row
    # can be easily computed by using a sliding window of
    # size of cisRange and counting NaNs there: 
    kernel_size = 2 * cis_range_size + 1
    count_NaN_kernel = np.ones(kernel_size, dtype='int')
    numnans_cisRange = np.convolve(
                            np.isnan(bias),
                            count_NaN_kernel,
                            mode='same')
    mask[numnans_cisRange > max_nan] = True

    # the mask is ready and is full compliant with the Tyler's
    # if,then,else implementation
    return mask


def main():

    parser=argparse.ArgumentParser(description='Get cis percent for window range of bin of Hi-C cooler file' +
    ' for each row of Hi-C matrix')
    parser.add_argument('-i', help='input cooler Hi-C file or cooler URI', type=str, required=True)
    parser.add_argument('-r', help='window range', type=int, default=6000000)
    parser.add_argument('-d', help='number of diagonals to ignore', type=int, default=0)
    parser.add_argument('-n', help='maximum percent of NAs allowed in window', type=float, default=10.0)
    args=parser.parse_args()

    # percent of NANs allowed
    p = args.n/100
    NaNsfraction = p

    # how many diagonals to ignore
    ignore_diags = args.d

    # genomic range - /2 to account for
    # upstream/downstream counting:
    grange = args.r/2

    # init a cooler object
    path_uri = args.i
    clr = cooler.Cooler(path_uri)

    # introducing binsize:
    binsize = clr.info['bin-size']
    range_bins=int(grange/binsize)

    # Number of NaN-row/cols we allow
    # a given CisRange band to intersect:
    max_NAs = 2 * range_bins * NaNsfraction

    # bins, number of bins, balancing weight, etc.
    n_bins = clr.info['nbins']
    clr_bins = clr.bins()[:]
    bias = clr_bins['weight'].values

    # we will create a mask to fill cisRange
    # with NaNs where mask is True:
    mask = _get_mask(n_bins,
                     max_NAs,
                     bias,
                     clr.chromnames,
                     clr.offset,
                     range_bins)
    # the mask is ready and is full compliant with the Tyler's
    # if,then,else implementation


    # start with actual calculations,
    # replacing "cp = mf.get_cis_percent_range(f, args.r, p)":

    # filling NaNs in the balancing weights with 0-s:
    _bias = np.nan_to_num(bias)
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

    # some final  modifications we need to do
    # to the cooler-pipeline generated result
    # to make it indistinguishable from Tyler's:
    cp = (cov-diag/2)
    cp[mask] = np.nan
    clr_bins['cisRange'] = cp

    # Write output
    # strip cool/mcool and uri part from input fname:
    cfname = path_uri.split("::")[0]
    cfname = '.'.join(cfname.split(".")[:-1])
    OUT_fname = cfname + '_range' + str(int(2*grange/1000000)) + 'Mb_cispercent.bedGraph'
    # Only using 22 autosomes and X ...
    clr_bins.to_csv(OUT_fname,
                    sep='\t',
                    na_rep='NA',
                    columns=["chrom","start","end",'cisRange'],
                    header=True,
                    index=False)

    # # something about excluding chrY chrM and other stuff
    # y_chrom_bin_start =  f['chr_bin_range'][:][23][0]
    # for i, b in enumerate(f['bin_positions'][:][:y_chrom_bin_start]):
    #     if b[0] == 22:
    #         chrom = 'chrX'
    #     else:
    #         chrom = 'chr' + str(b[0]+1)
    #     OUT.write(chrom+'\t'+str(b[1])+'\t'+str(b[2])+
    #     '\t'+ str(cp[i])+'\n')
    # OUT.close()


if __name__ == '__main__':
    main()
