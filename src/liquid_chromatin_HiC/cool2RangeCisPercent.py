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

def _get_chunk_range_cis_coverage(chunk, pixel_weight_key='count'):
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

def _get_chunk_coverage(chunk, pixel_weight_key='count'):
    '''
    Compute "coverage" for a given chunk of pixels.
    Used for calculating total sum of a Hi-C heatmap's row - Tyler style

    Parameters
    ----------
    chunk : dict of dict/pd.DataFrame
        A cooler chunk of pixels produced by the cooler split-apply-combine pipeline.
    pixel_weight_key: str
        The key of a pixel chunk to retrieve pixel weights.

    Returns
    -------
    cov : np.array n_bins
        A numpy array with the coverage.
    '''
    bins = chunk['bins']
    pixels = chunk['pixels']
    n_bins = len(bins['chrom'])
    cov = np.zeros(n_bins)
    pixel_weights = pixels[pixel_weight_key]

    cov += np.bincount(pixels['bin1_id'], weights=pixel_weights, minlength=n_bins)
    cov += np.bincount(pixels['bin2_id'], weights=pixel_weights, minlength=n_bins)

    # skip dividing by 2
    # and beware the main diagonal double counting
    return cov

def _get_chunk_diag(chunk, pixel_weight_key='count'):
    '''
    Extract diagonal values from a chunk of pixels.
    Used for adjusting cooler's coverage to Tyler's sums of rows.
    The two differ by a diagonal.

    Parameters
    ----------
    chunk : dict of dict/pd.DataFrame
        A cooler chunk of pixels produced by the cooler split-apply-combine pipeline.
    pixel_weight_key: str
        The key of a pixel chunk to retrieve pixel weights.

    Returns
    -------
    diag_vals : np.array n_bins
        A numpy array with the diagonal values.
    '''
    bins = chunk['bins']
    pixels = chunk['pixels']
    n_bins = len(bins['chrom'])
    diag_vals = np.zeros(n_bins)
    pixel_weights = pixels[pixel_weight_key]

    diag_mask = pixels['bin1_id'] == pixels['bin2_id']
    # only one summation
    diag_vals += np.bincount(pixels['bin1_id'], weights=pixel_weights*diag_mask, minlength=n_bins)

    return diag_vals



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
    parser.add_argument('-w', help='specify name of the column with balancing weights', type=str, default="weight")
    parser.add_argument('-t', help='divide cis-range by the total "coverage" (row-sum interpretation)', type=bool, default=False)
    args=parser.parse_args()

    # percent of NANs allowed
    p = args.n/100
    NaNsfraction = p

    # how many diagonals to ignore
    ignore_diags = args.d

    # divide by total to better match Tyler's:
    divide_by_total = args.t

    # weight column name to use in cooler:
    weight_column_name = args.w

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
    bias = clr_bins[weight_column_name].values

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
    cisrange_cov = (
        chunks
            .pipe(_zero_range, range_bins=range_bins)
            .pipe(_timesouterproduct, bias = _bias)
            .pipe(_get_chunk_range_cis_coverage, pixel_weight_key="balanced")
            .reduce( np.add, np.zeros(n_bins) )
    )
    cp = cisrange_cov

    # calculate the diagonal to adjust cooler's coverages to Tyler's sum of rows:
    if ignore_diags == 0:
        chunks = cooler.tools.split(clr, chunksize=int(1e7), map=map, use_lock=False)
        diag = (
            chunks
                .pipe(_timesouterproduct, bias = _bias)
                .pipe(_get_chunk_diag, pixel_weight_key="balanced")
                .reduce( np.add, np.zeros(n_bins) )
        )
        # afterwards covs-diag should be the same as clr.maxtrix()[:].sum(axis=0) ...
        cp = cp - diag

    # divide by total to match Tyler's results even closer
    # theoreticaly there is no need to divide to this in cooler,
    # because balancing implies sum of rows (+diag is not ignored) = 1.0
    # but due to the discrepancies between cooler/Tyler interpretation of coverage
    # we need this ...
    if divide_by_total:
        chunks = cooler.tools.split(clr, chunksize=int(1e7), map=map, use_lock=False)
        chunks = chunks.pipe(_zero_diags, n_diags=ignore_diags) if ignore_diags>0 else chunks
        total_cov = (
            chunks
                .pipe(_timesouterproduct, bias = _bias)
                .pipe(_get_chunk_coverage, pixel_weight_key="balanced")
                .reduce( np.add, np.zeros(n_bins) )
        )
        # turn total coverage to sum of row's:
        total_cov = total_cov if ignore_diags > 0 else total_cov - diag
        cp = cp / total_cov


    # apply masks and register result for output:
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


if __name__ == '__main__':
    main()
