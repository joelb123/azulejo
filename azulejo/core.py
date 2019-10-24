# -*- coding: utf-8 -*-
"""
Core logic for computing subtrees
"""
#
# standard library imports
#
import sys
from collections import Counter, OrderedDict
from datetime import datetime
from pathlib import Path
from itertools import chain, combinations
#
# third-party imports
#
import click
import networkx as nx
import pandas as pd
import numpy as np
from plumbum import local
#
# package imports
#
from . import cli, logger
from .common import *
#
# global constants
#
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
MAX_BINS = 10
DEFAULT_STEPS = 16

class ElapsedTimeReport:
    def __init__(self, name):
        self.start_time = datetime.now()
        self.name = name

    def elapsed(self, next_name):
        now = datetime.now()
        seconds = (now - self.start_time).total_seconds()
        report = '%s phase took %d seconds' %(self.name, seconds)
        self.start_time = now
        self.name = next_name
        return report


def read_synonyms(filepath):
    synonym_dict = {}
    try:
        df = pd.read_csv(filepath, sep='\t')
    except FileNotFoundError:
        print('Synonym tsv file "%s" does not exist' % substrpath)
    except pd.errors.EmptyDataError:
        print('Synonym tsv "%s" is empty' % substrpath)
    if len(df):
        if '#file' in df:
            df.drop('#file', axis=1, inplace=True)
        key = list(set(('Substr', 'Dups')
                       ).intersection(set(df.columns)))[0]
        for group in df.groupby('id'):
            synonym_dict[group[0]] = [substr for substr in group[1][key]]
    return synonym_dict


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



def get_fasta_ids(fasta):
    idset = set()
    with fasta.open() as f:
        for line in f:
            if line.startswith('>'):
                idset.add(line.split()[0][1:])
    return list(idset)


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


def parse_clusters(outdir,
                   identity,
                   delete=True,
                   count_clusters=True,
                   synonyms={}):
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
        if len(synonyms):
            syn_ids = set(ids).intersection(synonyms.keys())
            [ids.extend(synonyms[i]) for i in syn_ids]
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


def prettyprint_float(x, digits):
    format_string = '%.'+ '%d'%digits + 'f'
    return (format_string%x).rstrip('0').rstrip('.')


@cli.command()
@click.argument('seqfile')
@click.option('--identity','-i', default=0.,
              help='Minimum sequence identity (float, 0-1). [default: lowest]')
@click.option('--min_id_freq', '-m', default=0, show_default=True,
              help='Minimum frequency of ID components.')
@click.option('--delete/--no-delete', '-d/-n', is_flag=True, default=True,
              help='Delete primary output of clustering. [default: delete]')
@click.option('--write_ids/--no-write_ids', '-w', is_flag=True, default=False,
              help='Write file of ID-to-clusters. [default: delete]')
@click.option('--do_calc/--no-do_calc', '-c/-x', is_flag=True, default=True,
              help='Write file of ID-to-clusters. [default: do_calc]')
@click.option('--substrs', help='subpath to file of substrings. [default: none]')
@click.option('--dups', help='subpath to file of duplicates. [default: none]')
def usearch_cluster(seqfile,
                    identity,
                    delete=True,
                    write_ids=False,
                    do_calc=True,
                    min_id_freq=0,
                    substrs=None,
                    dups=None):
    """Cluster above a global sequence identity threshold"""
    global logger
    if identity == 1.0:
        digits = '10000'
    else:
        digits = ('%.4f'%identity)[2:]
    try:
        inpath, dirpath = get_paths_from_file(seqfile)
    except FileNotFoundError:
        logger.error('Input file "%s" does not exist!', seqfile)
        sys.exit(1)
    stem = inpath.stem
    dirpath = inpath.parent
    outname = stem + '-nr-%s' %digits
    outdir = '%s/'%outname
    logfile = '%s.log'%outname
    outfilepath = dirpath/outdir
    logfilepath = dirpath/logfile
    histfilepath = dirpath / ('%s-degreedist.tsv' %outname)
    gmlfilepath = dirpath / ('%s.gml' %outname)
    statfilepath = dirpath / ('%s-stats.tsv' %outname)
    anyfilepath = dirpath / ('%s-anyhist.tsv' %outname)
    allfilepath = dirpath / ('%s-allhist.tsv' %outname)
    idpath = dirpath/ ('%s-ids.tsv' %outname)
    logger.info('%s%% sequence identity output to %s*',
                prettyprint_float(identity*100, 2),
                outname)
    if not delete:
        logger.debug('Cluster files will be kept in %s and %s', logfile, outdir)
    if write_ids:
        logger.debug('File of cluster ID usage will be written to %s and %s',
                     anyfilepath, allfilepath)
    if not do_calc:
        if not logfilepath.exists():
            logger.error('Previous results must exists, rerun with --do_calc')
            sys.exit(1)
        logger.debug('Using previous results for calculation')
    if min_id_freq:
        logger.debug("Minimum number of times ID's must occur to be counted: %d",
                     min_id_freq)
    usearch = local['usearch']
    synonyms = {}
    if substrs is not None:
        logger.debug('using duplicates in %s', dirpath/dups)
        synonyms.update(read_synonyms(dirpath/substrs))
    if dups is not None:
        logger.debug('using duplicates in %s', dirpath/dups)
        synonyms.update(read_synonyms(dirpath/dups))
    timer = ElapsedTimeReport('usearch')
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
        with local.cwd(dirpath):
            calculate = usearch['-cluster_fast', seqfile,
                                '-id', identity,
                                '-clusters', outdir,
                                '-log', logfile]
            calculate()
            logger.debug(timer.elapsed('parse'))
    run_stat_dict = OrderedDict([('divergence', 1.-identity)])
    parse_usearch_log(logfilepath,
                      run_stat_dict)
    run_stats = pd.DataFrame(list(run_stat_dict.items()),
                             columns=['stat', 'val'])
    run_stats.set_index('stat', inplace=True)
    run_stats.to_csv(statfilepath, sep='\t')
    if delete:
        logfilepath.unlink()
    cluster_graph, clusters, ids, degrees, degree_counts, any_counts, all_counts =\
        parse_clusters(outfilepath,
                       identity,
                       delete=delete,
                       synonyms=synonyms)
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
    logger.debug(timer.elapsed('graph'))
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
    logger.debug(timer.elapsed('final'))
    return run_stat_dict, cluster_graph, cluster_hist, any_hist, all_hist


@cli.command()
@click.argument('seqfile')
@click.option('--steps','-s', default=DEFAULT_STEPS, show_default=True,
              help='# of steps from lowest to highest')
@click.option('--min_id_freq', '-m', default=0, show_default=True,
              help='Minimum frequency of ID components.')
@click.option('--substrs', help='subpath to file of substrings. [default: none]')
@click.option('--dups', help='subpath to file of duplicates. [default: none]')
def cluster_in_steps(seqfile,
                     steps,
                     min_id_freq=0,
                     substrs=None,
                     dups=None):
    """Cluster in steps from low to 100% sequence identity"""
    try:
        inpath, dirpath = get_paths_from_file(seqfile)
    except FileNotFoundError:
        logger.error('Input file "%s" does not exist!', seqfile)
        sys.exit(1)
    stat_path = dirpath / (inpath.stem + STATFILE_SUFFIX)
    any_path = dirpath / (inpath.stem + ANYFILE_SUFFIX)
    all_path = dirpath / (inpath.stem + ALLFILE_SUFFIX)
    logsteps = [1.] + list(1. - np.logspace(IDENT_LOG_MIN, IDENT_LOG_MAX, num=steps))
    logger.info('Clustering at %d levels from %s%% to %s%% global sequence identity',
                steps,
                prettyprint_float(min(logsteps)*100., 2),
                prettyprint_float(max(logsteps)*100., 2))
    stat_list = []
    all_frames = []
    any_frames = []
    for id_level in logsteps:
        stats, graph, hist, any, all = usearch_cluster.callback(seqfile,
                                                       id_level,
                                                       min_id_freq=min_id_freq,
                                                       substrs=substrs,
                                                       dups=dups)
        stat_list.append(stats)
        any_frames.append(any)
        all_frames.append(all)
    logger.info('Collating results on %s', seqfile)
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


@cli.command()
@click.argument('infile')
def clusters_to_histograms(infile):
    """Compute histograms from a tab-delimited cluster file"""
    try:
        inpath, dirpath = get_paths_from_file(infile)
    except FileNotFoundError:
        logger.error('Input file "%s" does not exist!', infile)
        sys.exit(1)
    histfilepath = dirpath/(inpath.stem + '-sizedist.tsv')
    clusters = pd.read_csv(dirpath/infile, sep='\t', index_col=0)
    cluster_counter = Counter()
    for cluster_id, group in clusters.groupby(['cluster']):
        #cluster_id = int(cluster_id.lstrip('cl'))
        cluster_counter.update({len(group): 1})
    logger.info('writing to %s', histfilepath)
    cluster_hist = pd.DataFrame(list(cluster_counter.items()),
                                columns=['size', 'clusts'])
    total_clusters = cluster_hist['clusts'].sum()
    cluster_hist['%clusts'] = cluster_hist['clusts'] * 100. / total_clusters
    cluster_hist['%genes'] = cluster_hist['clusts']*cluster_hist['size'] * 100. / len(clusters)
    cluster_hist.sort_values(['size'], inplace=True)
    cluster_hist.set_index('size', inplace=True)
    cluster_hist.to_csv(histfilepath, sep='\t', float_format='%06.3f')

