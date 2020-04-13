# -*- coding: utf-8 -*-
"""
data analysis and plotting
"""
import sys
from collections import OrderedDict
#
# third-party imports
#
import click
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
#
# package imports
#
from . import cli, logger
from .common import *
#
# Global constants
#
IDENT_LOG_MIN = -3
IDENT_LOG_MAX = 0
EPSILON = 0.000001
FILETYPE = 'png'
MAX_BINS = 10

def make_histogram(dist, name, log10=False, bins=None):
    # do histogram plot with kernel density estimate
    dist = dist[dist<10]
    mean = dist.mean()
    if log10:
        dist = np.log10(dist)
    #if len(dist) < MAX_BINS:
    #    bins = len(dist)
    #else:
    #    bins = MAX_BINS
    sns.distplot(dist,
                 bins=None,
                 rug=False,
                 kde=False,
                 norm_hist=True,
                 rug_kws={'color': 'b'},
                 kde_kws={'color': 'k',
                          'linewidth': 1,
                          'label': 'KDE'},
                 hist_kws={'histtype': 'step',
                           'linewidth': 2,
                           'alpha': 1,
                           'color': 'b',
                           #'range':(0,20)
                           }
                 )
    plt.title('%s histogram of %d values, mean=%.1f'
              % (name, len(dist), mean))
    if log10:
        plt.xlabel('log ' + name)
    else:
        plt.xlabel(name)
    plt.ylabel('Frequency')
    #for ext in PLOT_TYPES:
    #    plt.savefig(LOG_PATH / ('%s-histogram.' % (name.rstrip('%')) + ext),
    #                bbox_inches='tight')
    plt.show()
    #plt.close('all')


def tick_function(X):
    X = X*3.-3
    vals = [('%f'%v).rstrip('0').rstrip('.')
            for v in (1. - 10**X)*100.]
    ticks = ['%s%%'%v for v in vals]
    return ticks

def log_deriv(X,Y):
    logX = -1.0*np.log10(X+EPSILON)
    logY = np.log10(Y)
    return np.gradient(logY) / np.gradient(logX)

@cli.command()
@click.argument('instemlist')
def analyze_clusters(dirname,
                     instemlist,
                     label,
                     reference=None,
                     on_id=None,
                     match_type=None):
    if match_type is None:
        matches = ['all', 'any']
    else:
        matches = [match_type]
    uniques = {}
    divergence = {}
    dirpath = Path(dirname)
    div_dist = {'all': {'ref': 0.0},
                'any': {'ref': 0.0}}
    print('ref=',reference)
    for stem in instemlist:
        print('stem=', stem)
        paths = {'all': dirpath / (stem + ALLFILE_SUFFIX),
                 'any': dirpath / (stem + ANYFILE_SUFFIX),
                 'stat': dirpath / (stem + STATFILE_SUFFIX)}
        stats = pd.read_csv(paths['stat'], sep='\t', index_col=0)
        uniques[stem] = stats['unique_seqs'].iloc[0]
        divergence[stem] = stats['divergence']
        if on_id is None:
            div_dist['all'][stem] = log_deriv(divergence[stem],
                                              stats['clusters'])
            div_dist['any'][stem] = None
            if stem == reference:
                div_dist['all']['ref'] = div_dist['all'][stem]
                div_dist['any']['ref'] = None
        else:
            for match in ['any', 'all']:
                data = pd.read_csv(paths[match], sep='\t', index_col=0)
                try:
                    div_dist[match][stem] = log_deriv(divergence[stem],
                                             data.loc[on_id])
                except KeyError:# this label wasn't found
                    div_dist[match][stem] = None
                if stem == reference:
                    div_dist[match]['ref'] = div_dist[match][stem]
    #
    # Make the plots
    #
    plt.style.use('seaborn-whitegrid')
    axes = {}
    fig, ax = plt.subplots(len(matches),
                                   sharex=True)
    try:
        for axis, i in enumerate(ax):
            axes[matches[i]] = axis
            loweraxis = ax[1]
    except TypeError:
        axes[matches[0]] = ax
        loweraxis = ax
    for stem in instemlist:
        for match in matches:
            if div_dist[match][stem] is None:
                continue
            axes[match].plot(divergence[stem],
                     div_dist[match][stem] - div_dist[match]['ref'],
                     label='%s'%(stem.replace(label+'.','')))
                                #uniques[stem]/1000.))
    if reference is None:
        if on_id is None:
            title = '%s Divergence Distribution' %label
            outfilestem = '%s_divergence_dist.'%label
        else:
            title = '%s Divergence Distribution on "%s"' %(label, on_id)
            outfilestem = '%s_divergence_dist_%s.' %(label, on_id)
    else:
        if on_id is None:
            title = '%s_Differential Divergence Distribution vs. %s' %(label,
                                                                reference)
            outfilestem = '%s_divergence_dist_vs%s.' %(label, reference)
        else:
            title = '%s Differential Divergence Distribution on "%s" vs. %s'\
                 %(label, on_id, reference)
            outfilestem = '%s_divergence_dist_on_%s_vs_ref.' %(label,
                                                             on_id)
    if reference is None:
        fig.text(0.02, 0.5, 'Logarithmic Derivative on Clusters',
             ha='center',
             va='center',
             rotation='vertical')
    else:
        fig.text(0.02, 0.5, 'Logarithmic Derivative Difference on Clusters',
             ha='center',
             va='center',
             rotation='vertical')
    if len(matches) == 2:
        fig.text(0.5, 0.47,'All in Cluster',
                 ha='center', va='center')
        fig.text(0.5, 0.89, 'Any in Cluster',
                 ha='center', va='center')
    else:
        fig.text(0.5, 0.91,'%s in Cluster'%matches[0].capitalize(),
                 ha='center', va='center')
    loweraxis.set(xlabel='Divergence on Sequence Identity')
    loweraxis.legend(loc='upper left')
    fig.suptitle(title)
    plt.xscale('log')
    limits = [0.001,1.]
    new_tick_locations = np.array([0., 1./3., 2./3., 1.0])
    loweraxis.set_xlim(limits)
    axes['second'] = loweraxis.twiny()
    axes['second'].set_xlim(limits)
    axes['second'].set_xticks(new_tick_locations)
    axes['second'].set_xticklabels(tick_function(new_tick_locations))
    axes['second'].set_xlabel('   ')
    #r'%Identity')
    #plt.ylim([-0.002,0.002])
    outfilename = outfilestem + '%s'%FILETYPE
    print('saving plot to %s'%outfilename)
    plt.savefig(dirpath/outfilename, dpi=200)
    plt.show()

def do_cuts(obs, high, low, label):
    if high > 0.0:
        hicuts = obs[obs > high]
        obs = obs[obs <= high]
        if len(hicuts):
            hifilename = label + '_hicuts.tsv'
            logger.info('%d observations dropped by high-side cutoff of %.2f written to %s',
                        len(hicuts), high, hifilename)
            logger.info(hicuts)
    if low > 0.0:
        locuts = obs[obs > low]
        obs = obs[obs >= low]
        logger.info('%d observations dropped by low-side cutoff of %.2f',
                    len(locuts), low)
        if len(locuts):
            lofilename = label + '_locuts.tsv'
            logger.info('%d observations dropped by low-side cutoff of %.2f written to %s',
                        len(locuts), low, lofilename)
            logger.info(locuts)
    return obs

@cli.command()
@click.option('--hi_cutoff', default=2.0, show_default=True,
              help='Disregard above this value.')
@click.option('--lo_cutoff', default=0.0, show_default=True,
              help='Disregard below this value.')
@click.argument('cluster_size')
@click.argument('combinedfile')
def outlier_length_dist(hi_cutoff, lo_cutoff, cluster_size, combinedfile):
    """Plot the normalized length distribution of singletons in clusters"""
    cluster_size = int(cluster_size)
    if not cluster_size > 0:
        logger.error('Positive cluster size must be specified')
        sys.exit(1)
    clusters = pd.read_csv(combinedfile, sep='\t', index_col=0)
    norm_lengths = []
    for cluster_id, cluster in clusters.groupby('cluster'):
        if cluster['siz'].iloc[0] != cluster_size:
            # not the right size
            continue
        if len(set(cluster['sub'])) != 2:
            # not just two subclusters
            continue
        if 1 not in set(cluster['sub_siz']):
            # no singleton cluster
            continue
        singleton = cluster[cluster['sub_siz'] == 1]
        length = singleton['norm'].iloc[0]
        norm_lengths.append(length)
    norm_lengths = np.array(norm_lengths)
    norm_lengths = do_cuts(norm_lengths, hi_cutoff, lo_cutoff, 'len')
    logger.info('%d singleton outliers in clusters of size %d',
                len(norm_lengths), cluster_size)
    logger.info('min:\t%.3f', min(norm_lengths))
    logger.info('max:\t%.3f', max(norm_lengths))
    logger.info('mean: %.3f', norm_lengths.mean())
    ax = sns.distplot(norm_lengths, bins=100,
                      kde_kws={'label':'KDE'})
    ax.set_xlabel('Normalized Length of Singleton')
    plt.title('Length distribution of %d singleton subclusters'
              % (len(norm_lengths)))
    outfilename = 'norm_len_dist.%s'%FILETYPE
    logger.info('saving plot to %s', outfilename)
    #plt.yscale('log')
    plt.savefig(outfilename, dpi=200)
    plt.show()

@cli.command()
@click.option('--hi_cutoff', default=0.0, show_default=True,
              help='Disregard above this value.')
@click.option('--lo_cutoff', default=0.0, show_default=True,
              help='Disregard below this value.')
@click.argument('cluster_size')
@click.argument('combinedfile')
def length_std_dist(cluster_size, hi_cutoff, lo_cutoff, combinedfile):
    """Plot the normalized length distribution of singletons in clusters"""
    cluster_size = int(cluster_size)
    if not cluster_size > 0:
        logger.error('Positive cluster size must be specified')
        sys.exit(1)
    clusters = pd.read_csv(combinedfile, sep='\t', index_col=0)
    stds = []
    for cluster_id, cluster in clusters.groupby('cluster'):
        if cluster['siz'].iloc[0] != cluster_size:
            # not the right size
            continue
        if len(set(cluster['sub'])) != 1:
            # Only one subcluster
            continue
        val = cluster['std'].iloc[0]
        stds.append(val)
    stds = np.array(stds)
    pct_zeros = len(stds[stds == 0.0])*100/len(stds)
    stds = do_cuts(stds, hi_cutoff, lo_cutoff, 'stds')
    logger.info('%d single-subgroup clusters of size %d', len(stds), cluster_size)
    logger.info('%.1f %% zeroes, max is %.2f', pct_zeros, max(stds))
    logger.info('mean is %.3f', stds.mean())
    logbins = np.logspace(0.7, 3, 100)
    ax = sns.distplot(stds, bins=logbins, kde=False)
    ax.set_xlabel('Standard Deviation of Single-Subgroup Clusters')
    title = 'Length Standard Deviation distribution of %d clusters' %len(stds)
    plt.title(title)
    outfilename = 'std_dist.%s'%FILETYPE
    logger.info('saving plot to %s', outfilename)
    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(outfilename, dpi=200)
    plt.show()