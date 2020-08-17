"""
Saddle plot code.

Authors
~~~~~~~
* Anton Goloborodko
* Nezar Abdennur

"""
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.linalg import toeplitz
from cytoolz import merge
import numpy as np
import pandas as pd


pal = sns.color_palette('muted')


def digitize_track(
        bins,
        get_track,
        get_mask,
        chromosomes,
        prange,
        by_percentile=False):
    """
    Digitize per-chromosome genomic tracks.
    
    Parameters
    ----------
    get_track : function
        A function returning a genomic track given a chromosomal name/index.
    get_mask : function
        A function returning a binary mask of valid genomic bins given a 
        chromosomal name/index.
    chromosomes : list or iterator
        A list of names/indices of all chromosomes.
    bins : int or list or numpy.ndarray
        If list or numpy.ndarray, `bins` must contain the bin edges for 
        digitization.
        If int, then the specified number of bins will be generated 
        automatically from the genome-wide range of the track values.
    prange : pair of floats
        The percentile of the genome-wide range of the track values used to 
        generate bins. E.g., if `prange`=(2. 98) the lower bin would 
        start at the 2-nd percentile and the upper bin would end at the 98-th 
        percentile of the genome-wide signal.
        Use to prevent the extreme track values from exploding the bin range.
        Is ignored if `bins` is a list or a numpy.ndarray.
    by_percentile : bool
        If true then the automatically generated bins will contain an equal 
        number of genomic bins genome-wide (i.e. track values are binned 
        according to their percentile). Otherwise, bins edges are spaced equally. 
        Is ignored if `bins` is a list or a numpy.ndarray.
        
    Returns
    -------
    digitized : dict
        A dictionary of the digitized track, split by chromosome.
        The value of -1 corresponds to the masked genomic bins, the values of 0 
        and the number of bins correspond to the values lying below and above 
        the bin range limits. 
    binedges : numpy.ndarray
        The edges of bins used to digitize the track.
    """
    
    if not hasattr(bins, '__len__'):
        fulltrack = np.concatenate([
            get_track(chrom)[get_mask(chrom)] 
                for chrom in chromosomes
        ])
        
        if by_percentile:
            perc_edges = np.linspace(prange[0], prange[1], bins + 1)
            binedges = np.percentile(fulltrack, perc_edges)
        else:
            lo = np.percentile(fulltrack, prange[0])
            hi = np.percentile(fulltrack, prange[1])
            binedges = np.linspace(lo, hi, bins + 1)
    else:
        binedges = bins

    digitized = {}
    for chrom in chromosomes:
        x = np.digitize(get_track(chrom), binedges, right=False)        
        x[~get_mask(chrom)] = -1
        digitized[chrom] = x
        
    return digitized, binedges


def make_saddle(
        get_matrix,
        get_digitized,
        chromosomes,
        contact_type,
        verbose=False):
    """
    Make a matrix of average interaction probabilities between genomic bin pairs
    as a function of a specified genomic track. The provided genomic track must
    be pre-binned (i.e. digitized).
    
    Parameters
    ----------
    get_matrix : function
        A function returning an matrix of interaction between two chromosomes 
        given their names/indicies.
    get_digitized : function
        A function returning a track of the digitized target genomic track given
        a chromosomal name/index.
    chromosomes : list or iterator
        A list of names/indices of all chromosomes.    
    contact_type : str
        If 'cis' then only cis interactions are used to build the matrix.
        If 'trans', only trans interactions are used.
    verbose : bool
        If True then reports progress.

    Returns
    -------
    interaction_sum : 2D array
        The matrix of summed interaction probability between two genomic bins 
        given their values of the provided genomic track.
    interaction_count : 2D array
        The matrix of the number of genomic bin pairs that contributed to the 
        corresponding pixel of ``interaction_sum``.

    """
    if contact_type not in ['cis', 'trans']:
        raise ValueError("The allowed values for the contact_type "
                         "argument are 'cis' or 'trans'.")
    
    n_bins = max([
        get_digitized(chrom).max() 
            for chrom in chromosomes
    ])
    
    interaction_sum   = np.zeros((n_bins, n_bins))
    interaction_count = np.zeros((n_bins, n_bins))
    
    for k, chrom1 in enumerate(chromosomes):
        for chrom2 in chromosomes[k:]:
            if (((contact_type == 'trans') and (chrom1 == chrom2)) or 
                ((contact_type == 'cis') and (chrom1 != chrom2))):
                continue
                
            matrix = np.array(get_matrix(chrom1, chrom2))

            if verbose:
                print('chromosomes {} vs {}'.format(chrom1, chrom2))
                
            for i in range(n_bins):
                row_mask = (get_digitized(chrom1) == i)
                for j in range(n_bins):
                    col_mask = (get_digitized(chrom2) == j)
                    data = matrix[row_mask, :][:, col_mask]
                    data = data[np.isfinite(data)]
                    interaction_sum[i, j]   += np.sum(data)
                    interaction_count[i, j] += float(len(data))

    interaction_sum   += interaction_sum.T
    interaction_count += interaction_count.T
    
    return interaction_sum, interaction_count


def _saddleplot(binedges,
                digitized,
                saddledata,
                contact_type, 
                color,
                cbar_label=None,
                fig_kws=None,
                heatmap_kws=None, 
                margin_kws=None):
    """
    Plot saddle data and signal histograms in the margins.

    """
    n_bins = len(binedges) - 1
    lo, hi = 0, n_bins  #-0.5, n_bins - 1.5

    # Populate kwargs
    fig_kws = merge(
        dict(figsize=(5, 5)),
        fig_kws if fig_kws is not None else {}

    )

    heatmap_kws = merge(
        dict(aspect='auto', 
             cmap='coolwarm', 
             interpolation='none',
             vmin=-1,
             vmax=1),
        heatmap_kws if heatmap_kws is not None else {}, 
    )
    vmin = heatmap_kws['vmin']
    vmax = heatmap_kws['vmax']

    margin_kws = merge(
        dict(bins=n_bins,
             range=(0, len(binedges)),
             histtype='stepfilled',
             edgecolor='k',
             facecolor=color,
             linewidth=1),
        margin_kws if margin_kws is not None else {},
    )

    # layout
    gs = GridSpec(
        nrows=3, 
        ncols=3, 
        width_ratios=[0.2, 1, 0.1], 
        height_ratios=[0.2, 1, 0.1],
        wspace=0.05,
        hspace=0.05,
    )
    fig = plt.figure(**fig_kws)

    # heatmap
    ax = ax1 = plt.subplot(gs[4])
    img = ax.imshow(np.log10(saddledata), **heatmap_kws)
    plt.xticks(
        np.arange(0, n_bins)[::5],
        ['{:0.4f}'.format(t) 
            for t in ((binedges[1:] + binedges[:-1])/2)[::5]],
        rotation=90,
    )
    plt.yticks([])
    plt.xlim(lo, hi)
    plt.ylim(hi, lo)

    # left margin
    plt.subplot(gs[3])
    plt.hist(np.concatenate(list(digitized.values())), 
             **merge(margin_kws, {'orientation': 'horizontal'}))
    plt.xticks([])
    plt.yticks(
        np.arange(0, n_bins)[::5],
        ['{:0.4f}'.format(t) 
            for t in ((binedges[1:] + binedges[:-1])/2)[::5]],
    )
    plt.xlim(plt.xlim()[1], plt.xlim()[0])  # fliplr
    plt.ylim(hi, lo)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['bottom'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    # top margin
    plt.subplot(gs[1])
    plt.hist(np.concatenate(list(digitized.values())), **margin_kws)
    plt.xticks([])
    plt.yticks([])
    plt.xlim(lo, hi)
    plt.ylim(plt.ylim()[0], plt.ylim()[1])  # correct
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)

    # colorbar
    plt.subplot(gs[5])
    cb = plt.colorbar(
        img, 
        fraction=0.8, 
        label=cbar_label)
    if vmin is not None and vmax is not None:
        cb.set_ticks(np.arange(vmin, vmax + 0.001, 0.5))
    plt.grid(False)
    plt.axis('off')

    return fig



def make_cis_obsexp_fetcher(clr, expected, name):
    g = expected.groupby('chrom')
    return lambda chrom, _: (
            clr.matrix().fetch(chrom) / 
                toeplitz(g.get_group(chrom)[name].values)
        )


def make_trans_obsexp_fetcher(clr, trans_avg):
    return lambda chrom1, chrom2: (
            clr.matrix().fetch(chrom1, chrom2) / trans_avg
        )


def make_track_fetcher(df, name):
    g = df.groupby('chrom')
    return lambda chrom: g.get_group(chrom)[name]


def saddleplot(
        clr,
        track,
        track_name,
        contact_type='cis',
        expected=None,
        expected_name=None,
        chromosomes=None,
        n_bins=50,
        prange=(0.25, 99.75),
        by_percentile=False,
        color=pal[2],
        fig_kws=None,
        heatmap_kws=None, 
        margin_kws=None):
    """
    Make a saddle plot for an arbitrary signal track on the genomic bins of 
    a contact matrix.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object
    track : pandas.DataFrame
        bedGraph-like data frame ('chrom', 'start', 'end', <name>).
    track_name : str
        Name of value column in ``track``.
    n_bins : int, optional
        Number of bins for digitizing the signal track.
    chromosomes : list of str, optional
        Restricted set of chromosomes to use. Default is to use all chromosomes
        in the ``track`` dataframe.
    contact_type : 'cis' or 'trans'
        Type of contacts to aggregate. Default is cis.
    expected : pandas.DataFrame or number, optional
        Expected contact frequency. If ``contact_type`` is 'cis', then this is
        a data frame of ('chrom', 'distance', <name>). If ``contact_type`` is 
        'trans', then this is a single precomputed number equal to the average 
        trans contact frequency.
    expected_name : str, optional
        Required if ``contact_type`` is 'cis'. Name of value column in the
        ``expected`` data frame.
    by_percentile : bool, optional
        Whether to bin the signal track by percentile rather than by value.
        Default is False.
    bin_range_percentile : float, optional
        ...
    heatmap_kws : dict, optional
        Extra keywords to pass to ``imshow`` for saddle heatmap.
    margin_kws : dict, optional
        Extra keywords to pass to ``hist`` for left and top margins.

    Returns
    -------
    figure

    """
    if contact_type == 'cis':
        if not isinstance(expected, pd.DataFrame):
            raise ValueError("For cis contacts, `expected` must be a "
                             "bedGraph-like dataframe")
        if expected_name is None:
            raise ValueError("Need a value column name for `expected`")

        obsexp_func = make_cis_obsexp_fetcher(clr, expected, expected_name)

    elif contact_type == 'trans':
        if not np.isscalar(expected):
            raise ValueError("For trans contacts, `expected` must be the "
                             "scalar average trans value")

        obsexp_func = make_trans_obsexp_fetcher(clr, expected)
    else:
        raise ValueError("`contact_type` must be 'cis' or 'trans'")

    if chromosomes is None:
        chromosomes = list(track['chrom'].unique())

    track_fetcher = make_track_fetcher(track, track_name)
    track_mask_fetcher = lambda chrom: (
        (track.groupby('chrom')
              .get_group(chrom)[track_name]
              .notnull()) &
        (track.groupby('chrom')
              .get_group(chrom)[track_name] != 0)
    ).values

    # Bin the track
    digitized, binedges = digitize_track(
        n_bins,
        track_fetcher,
        track_mask_fetcher,
        chromosomes,
        prange=prange,
        by_percentile=by_percentile
    )

    # Aggregate contacts
    sum_, count = make_saddle(
        obsexp_func, 
        lambda chrom: digitized[chrom],
        chromosomes,
        contact_type,
        verbose=False)
    saddledata = sum_ / count

    # Plotting
    fig = _saddleplot(
        binedges,
        digitized,
        saddledata,
        contact_type, 
        color,
        cbar_label='log10 (contact frequency / mean) in {}'.format(contact_type),
        fig_kws=None,
        heatmap_kws=None, 
        margin_kws=None)

    return binedges, digitized, saddledata, fig