# -*- coding: utf-8 -*-
'''
azulejo -- tile phylogenetic space with subtrees
'''
import functools
import locale
import logging
import os
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path
#
# third-party imports
#
import click
import coverage
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
#
# package imports
#
from .version import version as VERSION
from .common import *
#
# start coverage
#
coverage.process_startup()
#
# set locale so grouping works
#
for localename in ['en_US', 'en_US.utf8', 'English_United_States']:
    try:
        locale.setlocale(locale.LC_ALL, localename)
        break
    except locale.Error:
        continue
#
# global constants
#
PROGRAM_NAME='azulejo'
AUTHOR = 'Joel Berendzen'
EMAIL = 'joelb@ncgr.org'
COPYRIGHT = """Copyright (C) 2019, NCGR. All rights reserved.
"""
LOG_DIR = 'logs'
LOG_PATH = Path('.') / LOG_DIR
# defaults for command line
DEFAULT_FILE_LOGLEVEL = logging.DEBUG
DEFAULT_STDERR_LOGLEVEL = logging.INFO

IDENT_LOG_MIN = -3
IDENT_LOG_MAX = 0
EPSILON = 0.000001
FILETYPE = 'pdf'
MAX_BINS = 10
#
# Class definitions.
#
class CleanInfoFormatter(logging.Formatter):
    def __init__(self, fmt='%(levelname)s: %(message)s'):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        if record.levelno == logging.INFO:
            return record.getMessage()
        return logging.Formatter.format(self, record)
#
# time stamp for start
#
STARTTIME = datetime.now()
#
# global logger object
#
logger = logging.getLogger(PROGRAM_NAME)
#
# private context function
#
_ctx = click.get_current_context
#
# Helper functions
#
def init_dual_logger(file_log_level=DEFAULT_FILE_LOGLEVEL,
                     stderr_log_level=DEFAULT_STDERR_LOGLEVEL):
    '''Log to stderr and to a log file at different levels
    '''
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            global logger
            # find out the verbose/quiet level
            if _ctx().params['verbose']:
                _log_level = logging.DEBUG
            elif _ctx().params['quiet']:
                _log_level = logging.ERROR
            else:
                _log_level = stderr_log_level
            logger.setLevel(file_log_level)
            stderrHandler = logging.StreamHandler(sys.stderr)
            stderrFormatter = CleanInfoFormatter()
            stderrHandler.setFormatter(stderrFormatter)
            stderrHandler.setLevel(_log_level)
            logger.addHandler(stderrHandler)
            if _ctx().params['log']:  # start a log file in LOG_PATH
                logfile_path = LOG_PATH / (PROGRAM_NAME + '.log')
                if not LOG_PATH.is_dir():  # create LOG_PATH
                    try:
                        logfile_path.parent.mkdir(mode=0o755, parents=True)
                    except OSError:
                        logger.error('Unable to create log directory "%s"',
                                     logfile_path.parent)
                        raise OSError
                else:
                    if logfile_path.exists():
                        try:
                            logfile_path.unlink()
                        except OSError:
                            logger.error('Unable to remove log file "%s"',
                                         logfile_path)
                            raise OSError
                logfileHandler = logging.FileHandler(str(logfile_path))
                logfileFormatter = logging.Formatter(
                    '%(levelname)s: %(message)s')
                logfileHandler.setFormatter(logfileFormatter)
                logfileHandler.setLevel(file_log_level)
                logger.addHandler(logfileHandler)
            logger.debug('Command line: "%s"', ' '.join(sys.argv))
            logger.debug('%s version %s', PROGRAM_NAME, VERSION)
            logger.debug('Run started at %s', str(STARTTIME)[:-7])
            return f(*args, **kwargs)
        return wrapper
    return decorator


def init_user_context_obj(initial_obj=None):
    '''Put info from global options into user context dictionary
    '''
    def decorator(f):
        @functools.wraps(f)
        def wrapper(*args, **kwargs):
            global config_obj
            if initial_obj is None:
                _ctx().obj = {}
            else:
                _ctx().obj = initial_obj
            ctx_dict = _ctx().obj
            if _ctx().params['verbose']:
                ctx_dict['logLevel'] = 'verbose'
            elif _ctx().params['quiet']:
                ctx_dict['logLevel'] = 'quiet'
            else:
                ctx_dict['logLevel'] = 'default'
            return f(*args, **kwargs)
        return wrapper
    return decorator

@click.group(epilog=AUTHOR + ' <' + EMAIL + '>. ' + COPYRIGHT)
@click.option('--warnings_as_errors', is_flag=True, show_default=True,
              default=False, help='Warnings cause exceptions.')
@click.option('-v', '--verbose', is_flag=True, show_default=True,
              default=False, help='Debug info to stderr.')
@click.option('-q', '--quiet', is_flag=True, show_default=True,
              default=False, help='Suppress logging to stderr.')
@click.option('--log/--no-log', is_flag=True, show_default=True,
              default=True, help='Write analysis in ./' + LOG_DIR + '.')
@click.version_option(version=VERSION, prog_name=PROGRAM_NAME)
@init_dual_logger()
@init_user_context_obj()
def cli(warnings_as_errors, verbose, quiet, log):
    """azulejo -- tiling gene in subtrees across phylogenetic space

    If COMMAND is present, and --no_log was not invoked,
    a log file named azulejo-COMMAND.log
    will be written in the ./logs/ directory.
    """
    if warnings_as_errors:
        logger.debug('Runtime warnings (e.g., from pandas) will cause exceptions')
        warnings.filterwarnings('error')

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


def cluster_in_steps(infile,
                     steps,
                     min_id_freq=0,
                     substrings=None,
                     duplicates=None):
    inpath = Path(infile)
    dirpath = inpath.parent
    stat_path = dirpath / (inpath.stem + STATFILE_SUFFIX)
    any_path = dirpath / (inpath.stem + ANYFILE_SUFFIX)
    all_path = dirpath / (inpath.stem + ALLFILE_SUFFIX)
    logsteps = [1.] + list(1. - np.logspace(IDENT_LOG_MIN, IDENT_LOG_MAX, num=steps))
    print('clustering %s at %d levels from %f to %f sequence identity'
          %(infile, steps, min(logsteps), max(logsteps)))
    stat_list = []
    all_frames = []
    any_frames = []
    for id_level in logsteps:
        stats, graph, hist, any, all = usearch_cluster(infile,
                                                       id_level,
                                                       min_id_freq=min_id_freq,
                                                       substrings=substrings,
                                                       duplicates=duplicates)
        stat_list.append(stats)
        any_frames.append(any)
        all_frames.append(all)
    print('collating results')
    #
    # Concatenate and write stats
    #
    stats = pd.DataFrame(stat_list)
    stats.to_csv(stat_path, sep='\t')
    #
    # Concatenate any/all data
    #
    any = pd.concat(any_frames, axis=1, join='inner',
                    sort=True, ignore_index=False)
    any.to_csv(any_path, sep='\t')
    all = pd.concat(all_frames, axis=1, join='inner',
                    sort=True, ignore_index=False)
    all.to_csv(all_path, sep='\t')


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
    for stem in instemlist:
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
            outfilestem = '%s_divergence_dist_on_%s_vs_%s.' %(label,
                                                             on_id,
                                                             reference)
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

def clusters_to_histograms(dirname, infile):
    dirpath = Path(dirname)
    infilepath = Path(infile)
    histfilepath = dirpath/(infilepath.stem + '-sizedist.tsv')
    clusters = pd.read_csv(dirpath/infile, sep='\t', index_col=0)
    cluster_counter = Counter()
    for cluster_id, group in clusters.groupby(['cluster']):
        cluster_id = int(cluster_id.lstrip('cl'))
        cluster_counter.update({len(group): 1})
    print('writing to %s'%histfilepath)
    cluster_hist = pd.DataFrame(list(cluster_counter.items()),
                                columns=['size', 'clusts'])
    total_clusters = cluster_hist['clusts'].sum()
    cluster_hist['%clusts'] = cluster_hist['clusts'] * 100. / total_clusters
    cluster_hist['%genes'] = cluster_hist['clusts']*cluster_hist['size'] * 100. / len(clusters)
    cluster_hist.sort_values(['size'], inplace=True)
    cluster_hist.set_index('size', inplace=True)
    cluster_hist.to_csv(histfilepath, sep='\t', float_format='%06.3f')


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

# core logic imports here
from .core import usearch_cluster

def test():
    TESTDIR = '/home/localhost/joelb/preclust/'
    genomes = {'medtr': {'all': ['jemalong_A17', 'R108_HM340', 'HM004',
                                 'HM010', 'HM022', 'HM023',
                                 'HM034', 'HM050', 'HM056',
                                 'HM058', 'HM060', 'HM095',
                                 'HM125', 'HM129', 'HM185',
                                 'HM324',
                                 'jema+HM056',
                                 'jema+R108',
                                 'all'
                                ],
                         'refs': ['jemalong_A17']},
               'glyma':{'all': ['Williams',
                                'Lee',
                                'Zhonghuang',
                                'all',
                                #'Wm+glyso',
                                #'Wm+phavu',
                                #'Wm+jemalong',
                                ],
                        'ref': 'Williams'},
               'glyso': {'all': ['W05',
                                 'PI'],
                         'ref': 'W05'},
               'glyma+glyso': {'all': ['all'],
                               'ref': 'all'}
               }
    species = 'glyma+glyso'
    #operation = 'compute'
    operation = 'single'
    #operation = 'plot'
    #operation = 'condense'
    #operation = 'compare'
    print('doing species %s' %species)
    if operation == 'compute':
        for gnm in genomes[species]['all']:
            cluster_in_steps(TESTDIR + species + '/' + gnm,
                             '%s.%s.faa'%(species, gnm),
                             16,
                             substrs='log/Substr.tsv',
                             dups='log/Dups.tsv',
                             min_id_freq=10)
    elif operation == 'single':
        gnm = genomes[species]['all'][0]
        usearch_cluster(TESTDIR + species + '/' +gnm +'/%s.%s.faa'%(species, gnm),
                        0.984,
                        delete=False,
                        write_ids=True,
                        #do_calc=False,
                        subtrs='log/Substr.tsv',
                        dupss='log/Dups.tsv',
                        min_id_freq=10,
                        )
    elif operation == 'plot':
        namelist = []
        for gnm in genomes[species]['all']:
            namelist.append('%s.%s'%(species,gnm))
        analyze_clusters(TESTDIR + species,
                        namelist, species,
                        reference='%s.%s'%(species,
                                          genomes[species]['ref']),
                         on_id='glyma',
                         match_type='all')
    elif operation == 'condense':
        clusters_to_histograms(TESTDIR+'glyma+glyso/steven',
                               'steven_clusters.tsv')
    elif operation == 'compare':
        compare_clusters(TESTDIR+'glyma+glyso/steven/steven_clusters.tsv',
                         TESTDIR+'glyma+glyso/all/glyma+glyso.all-nr-984-ids.tsv')
