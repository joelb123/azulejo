# -*- coding: utf-8 -*-
"""
data analysis and plotting
"""
import os
import sys
from collections import Counter
from pathlib import Path
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

#@cli.command()
#@click.argument('instemlist')
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


def compare_clusters(file1, file2):
    path1 = Path(file1)
    path2 = Path(file2)
    commondir = Path(os.path.commonpath([path1, path2]))
    missing1 = commondir/'notin1.tsv'
    missing2 = commondir/'notin2.tsv'
    clusters1 = pd.read_csv(path1, sep='\t', index_col=0)
    print('%d members in %s'%(len(clusters1), file1))
    clusters2 = pd.read_csv(path2, sep='\t', index_col=0)
    print('%d members in %s'%(len(clusters2), file2))
    ids1 = set(clusters1['id'])
    ids2 = set(clusters2['id'])
    notin1 = pd.DataFrame(ids2.difference(ids1), columns=['id'])
    notin1.sort_values('id', inplace=True)
    notin1.to_csv(missing1, sep='\t')
    notin2 = pd.DataFrame(ids1.difference(ids2), columns=['id'])
    notin2.sort_values('id', inplace=True)
    notin2.to_csv(missing2, sep='\t')

    print('%d ids not in ids1' %len(notin1))
    print('%d ids not in ids2' %len(notin2))
    print('%d in %s after dropping'%(len(clusters1), file1))
    #print(notin2)
