# -*- coding: utf-8 -*-
'''
azulejo -- tile phylogenetic space with subtrees
'''
import logging
import zlib
from collections import Counter, OrderedDict
from operator import itemgetter
from pathlib import Path
from itertools import chain, combinations
# third-party imports
from plumbum import local
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

UNITS = {'Mb':{'factor': 1,
               'outunits': 'MB'},
         'Gb':{'factor': 1024,
               'outunits': 'MB'},
         's':{'factor': 1,
              'outunits': 's'},
         'm':{'factor': 60,
              'outunits': 's'},
         'h':{'factor': 3600,
              'outunits': 's'}}
SEQ_IN_LINE = 6
IDENT_STATS_LINE = 7
FIRST_LOG_LINE = 14
LAST_LOG_LINE = 23
STAT_SUFFIXES = ['size', 'mem', 'time', 'memory']
RENAME_STATS = {'throughput':'throughput_seq_s',
                'time': 'CPU_time',
                'max_size': 'max_cluster_size',
                'avg_size': 'avg_cluster_size',
                'min_size': 'min_cluster_size',
                'seqs': 'unique_seqs',
                'singletons': 'singleton_clusters'}
ID_SEPARATOR = '.'
IDENT_LOG_MIN = -3
IDENT_LOG_MAX = 0
EPSILON = 0.000001
FILETYPE = 'pdf'
DEFAULT_ID_FREQ_CUTOFF = 10
__NAME__ = 'azulejo'
STATFILE_SUFFIX = '-%s_stats.tsv' %__NAME__
ANYFILE_SUFFIX = '-%s_ids-any.tsv' %__NAME__
ALLFILE_SUFFIX = '-%s_ids-all.tsv' %__NAME__
CLUSTFILE_SUFFIX = '-%s_clusts.tsv' %__NAME__
MAX_BINS = 10

def parse_usearch_log(filepath, rundict):
    with filepath.open() as logfile:
        for lineno, line in enumerate(logfile):
            if lineno < FIRST_LOG_LINE:
                if lineno == SEQ_IN_LINE:
                    split = line.split()
                    rundict['seqs_in'] = int(split[0])
                    rundict['singleton_seqs_in'] = int(split[4])
                if lineno == IDENT_STATS_LINE:
                    split = line.split()
                    rundict['max_identical_seqs'] = int(split[6].rstrip(','))
                    rundict['avg_identical_seqs'] = float(split[8])
                continue
            if lineno > LAST_LOG_LINE:
                break
            split = line.split()
            if split:
                stat = split[0].lower()
                if split[1] in STAT_SUFFIXES:
                    stat += '_' + split[1]
                    val = split[2]
                else:
                    val = split[1].rstrip(',')
                # rename poorly-named stats
                if stat in RENAME_STATS:
                    stat = RENAME_STATS[stat]
                # strip stats with units at the end
                conversion_factor = 1
                for unit in UNITS.keys():
                    if val.endswith(unit):
                        val = val.rstrip(unit)
                        conversion_factor = UNITS[unit]['factor']
                        stat += '_' + UNITS[unit]['outunits']
                        break
                # convert string values to int or float where possible
                try:
                    val = int(val)
                    val *= conversion_factor
                except:
                    try:
                        val = float(val)
                        val *= conversion_factor
                    except:
                        pass
                rundict[stat] = val
            #print(lineno, split)

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


def get_fasta_ids(fasta):
    idset = set()
    with fasta.open() as f:
        for line in f:
            if line.startswith('>'):
                idset.add(line.split()[0][1:])
    return idset

def parse_chromosome(id):
    #
    # If id contains an underscore, work on the
    # last part only (e.g., MtrunA17_Chr4g0009691)
    #
    undersplit = id.split('_')
    if len(undersplit) > 1:
        id = undersplit[-1].upper()
        if id.startswith('CHR'):
            id = id[3:]
    #
    # Chromosome numbers are integers suffixed by 'G'
    #
    try:
        chr = 'Chr'+str(int(id[:id.index('G')]))
    except:
        chr = None
    return chr

def parse_subids(id):
    subids = id.split(ID_SEPARATOR)
    subids += [chr for chr in [parse_chromosome(id) for id in subids]
                   if chr is not None]
    return subids

def parse_clusters(outdir, identity, delete=True, count_clusters=True):
    cluster_list = []
    id_list = []
    degree_list = []
    degree_counter = Counter()
    any_counter = Counter()
    all_counter = Counter()
    graph = nx.Graph()
    for fasta in outdir.glob('*'):
        cluster_id = int(fasta.name)

        ids = get_fasta_ids(fasta)
        n_ids = len(ids)
        degree_list.append(n_ids)
        degree_counter.update({n_ids:1})
        id_list += ids
        cluster_list += [cluster_id] * n_ids
        #
        # Do 'any' and 'all' counters
        #
        id_counter = Counter()
        id_counter.update(chain.from_iterable(
            [parse_subids(id) for id in ids]))
        if count_clusters:
            any_counter.update(id_counter.keys())
            all_counter.update([id for id in id_counter.keys()
                                if id_counter[id] == n_ids])
        elif n_ids > 1:
            any_counter.update({s: n_ids for s in id_counter.keys()})
            all_counter.update({id: n_ids for id in id_counter.keys()
                                if id_counter[id] == n_ids})
        #
        # Do graph components
        #
        if n_ids > 1:
            graph.add_nodes_from(ids)
            edges = combinations(ids, 2)
            graph.add_edges_from(edges, weight=identity)
        if delete:
            fasta.unlink()
    if delete:
        outdir.rmdir()
    return graph, cluster_list, id_list, degree_list, degree_counter, any_counter, all_counter


def usearch_cluster(dir,
                    seqfile,
                    identity,
                    delete=True,
                    write_ids=False,
                    do_calc=True,
                    min_id_freq=0):
    """cluster sequences using usearch, collecting stats"""
    print('   clustering %s by usearch at %f sequence identity'
          %(seqfile, identity))
    usearch = local['usearch']
    if identity == 1.0:
        digits = '1000'
    else:
        digits = ('%.3f'%identity)[2:]
    stem = Path(seqfile).stem
    outname = stem + '-nr-%s' %digits
    outdir = '%s/'%outname
    logfile = '%s.log'%outname
    dirpath = Path(dir)
    outfilepath = dirpath/outdir
    logfilepath = dirpath/logfile
    histfilepath = dirpath / ('%s-degreedist.tsv' %outname)
    gmlfilepath = dirpath / ('%s.gml' %outname)
    statfilepath = dirpath / ('%s-stats.tsv' %outname)
    anyfilepath = dirpath / ('%s-anyhist.tsv' %outname)
    allfilepath = dirpath / ('%s-allhist.tsv' %outname)
    idpath = dirpath/ ('%s-ids.tsv' %outname)
    if do_calc:
        #
        # Delete previous results, if any.
        #
        if outfilepath.exists() and outfilepath.is_file():
            outfilepath.unlink()
        elif outfilepath.exists() and outfilepath.is_dir():
            for file in outfilepath.glob('*'):
                file.unlink()
                pass
        else:
            outfilepath.mkdir()
        #
        # Do the calculation.
        #
        with local.cwd(dir):
            calculate = usearch['-cluster_fast', seqfile,
                                '-id', identity,
                                '-clusters', outdir,
                                '-log', logfile]
            calculate()
    run_stat_dict = OrderedDict([('divergence', 1.-identity)])
    parse_usearch_log(logfilepath, run_stat_dict)
    run_stats = pd.DataFrame(list(run_stat_dict.items()),
                             columns=['stat', 'val'])
    run_stats.set_index('stat', inplace=True)
    run_stats.to_csv(statfilepath, sep='\t')
    if delete:
        logfilepath.unlink()
    cluster_graph, clusters, ids, degrees, degree_counts, any_counts, all_counts =\
        parse_clusters(outfilepath, identity, delete=delete)
    #
    # Write out list of clusters and ids.
    #
    id_frame = pd.DataFrame.from_dict({'id': ids,
                            'cluster': clusters})
    id_frame.sort_values('cluster', inplace=True)
    id_frame = id_frame.reindex(['cluster', 'id'], axis=1)
    id_frame.reset_index(inplace=True)
    id_frame.drop(['index'], axis=1, inplace=True)
    id_frame.to_csv(idpath, sep='\t')
    del ids, clusters, id_frame
    #
    # Write out degree distribution.
    #
    cluster_hist = pd.DataFrame(list(degree_counts.items()), columns=['degree', 'clusters'])
    #cluster_hist['degree'] = cluster_hist['degree'] - 1
    cluster_hist.sort_values(['degree'], inplace=True)
    cluster_hist.set_index('degree', inplace=True)
    total_clusters = cluster_hist['clusters'].sum()
    cluster_hist['pct_total'] =   cluster_hist['clusters']*100./total_clusters
    cluster_hist.to_csv(histfilepath, sep='\t', float_format='%06.3f')
    del degree_counts
    #
    # Do histograms of "any" and "all" id usage in cluster
    #
    hist_value = '%f'%identity
    any_hist = pd.DataFrame(list(any_counts.items()), columns=['id', hist_value])
    any_hist.set_index('id', inplace=True)
    any_hist.sort_values(hist_value, inplace=True, ascending=False)
    all_hist = pd.DataFrame(list(all_counts.items()), columns=['id', hist_value])
    all_hist.set_index('id', inplace=True)
    all_hist.sort_values(hist_value, inplace=True, ascending=False)
    if min_id_freq:
        any_hist = any_hist[any_hist[hist_value] > min_id_freq]
        all_hist = all_hist[all_hist[hist_value] > min_id_freq]
    if write_ids:
        any_hist.to_csv(anyfilepath, sep='\t')
        all_hist.to_csv(allfilepath, sep='\t')
    #
    # Compute cluster stats
    #
    #degree_sequence = sorted([d for n, d in cluster_graph.degree()], reverse=True)
    #degreeCount = Counter(degree_sequence)
    #degree_hist = pd.DataFrame(list(degreeCount.items()),
    #                           columns=['degree', 'count'])
    #degree_hist.set_index('degree', inplace=True)
    #degree_hist.sort_values('degree', inplace=True)
    #degree_hist.to_csv(histfilepath, sep='\t')
    nx.write_gml(cluster_graph, gmlfilepath)
    return run_stat_dict, cluster_graph, cluster_hist, any_hist, all_hist


def cluster_in_steps(dirname, infile, steps, min_id_freq=0):
    stat_path = Path(dirname) / (Path(infile).stem + STATFILE_SUFFIX)
    any_path = Path(dirname) / (Path(infile).stem + ANYFILE_SUFFIX)
    all_path = Path(dirname) / (Path(infile).stem + ALLFILE_SUFFIX)
    logsteps = [1.] + list(1. - np.logspace(IDENT_LOG_MIN, IDENT_LOG_MAX, num=steps))
    print('clustering %s at %d levels from %f to %f sequence identity'
          %(infile, steps, min(logsteps), max(logsteps)))
    stat_list = []
    all_frames = []
    any_frames = []
    for id_level in logsteps:
        stats, graph, hist, any, all = usearch_cluster(dirname,
                                                       infile,
                                                       id_level,
                                                       min_id_freq=min_id_freq)
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
    histfilepath = dirpath/(infilepath.name + '-degreedist.tsv')
    clusters = pd.read_csv(dirpath/infile, sep='\t', index_col=0)

    cluster_counter = Counter()
    for cluster_id, group in clusters.groupby(['cluster']):
        cluster_id = int(cluster_id.lstrip('cl'))
        cluster_counter.update({len(group): 1})
    print('writing to %s'%histfilepath)
    cluster_hist = pd.DataFrame(list(cluster_counter.items()),
                                columns=['degree', 'clusters'])
    cluster_hist.sort_values(['degree'], inplace=True)
    cluster_hist.set_index('degree', inplace=True)
    total_clusters = cluster_hist['clusters'].sum()
    cluster_hist['pct_total'] = cluster_hist['clusters'] * 100. / total_clusters
    cluster_hist.to_csv(histfilepath, sep='\t', float_format='%06.3f')


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
                                #'Lee',
                                #'Zhonghuang',
                                'Wm+Lee',
                                'Wm+Zh',
                                #'all+glyso',
                                #'all',
                                #'Wm+glyso',
                                'Wm+glysoW05',
                                'Wm+glysoPI'
                                #'Wm+phavu',
                                #'Wm+jemalong',
                                ],
                        'ref': 'Williams'}}
    species = 'glyma'
    #operation = 'compute'
    operation = 'single'
    #operation = 'plot'
    #operation = 'condense'
    print('doing species %s' %species)
    if operation == 'compute':
        for gnm in genomes[species]['all']:
            cluster_in_steps(TESTDIR + species,
                    '%s.%s.faa'%(species, gnm), 16,
                             min_id_freq=DEFAULT_ID_FREQ_CUTOFF)
    elif operation == 'single':
        #gnm = genomes[species]['all'][0]
        gnm = 'all+glyso'
        usearch_cluster(TESTDIR + species,
                        '%s.%s.faa'%(species, gnm),
                        0.96,
                        delete=False,
                        write_ids=True,
                        #do_calc=False,
                        min_id_freq=DEFAULT_ID_FREQ_CUTOFF)
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
        clusters_to_histograms(TESTDIR+'glyma/steven',
                               'steven_clusters.tsv')
