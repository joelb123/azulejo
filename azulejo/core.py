# -*- coding: utf-8 -*-
"""
Core logic for computing subtrees
"""

# standard library imports
import contextlib
import os
import sys
from collections import Counter
from collections import OrderedDict
from datetime import datetime
from itertools import chain
from itertools import combinations

# third-party imports
import click
import dask.bag as db
import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import IUPACData
from dask.diagnostics import ProgressBar

# first-party imports
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import *

#
# global constants
#
UNITS = {
    "Mb": {"factor": 1, "outunits": "MB"},
    "Gb": {"factor": 1024, "outunits": "MB"},
    "s": {"factor": 1, "outunits": "s"},
    "m": {"factor": 60, "outunits": "s"},
    "h": {"factor": 3600, "outunits": "s"},
}
SEQ_IN_LINE = 6
IDENT_STATS_LINE = 7
FIRST_LOG_LINE = 14
LAST_LOG_LINE = 23
STAT_SUFFIXES = ["siz", "mem", "time", "memory"]
RENAME_STATS = {
    "throughput": "throughput_seq_s",
    "time": "CPU_time",
    "max_size": "max_cluster_size",
    "avg_size": "avg_cluster_size",
    "min_size": "min_cluster_size",
    "seqs": "unique_seqs",
    "singletons": "singleton_clusters",
}
ID_SEPARATOR = "."
IDENT_LOG_MIN = -3
IDENT_LOG_MAX = 0
EPSILON = 0.000001
FILETYPE = "pdf"
MAX_BINS = 10
DEFAULT_STEPS = 16
ALPHABET = IUPACData.protein_letters + "X" + "-"
#
# Classes
#
class ElapsedTimeReport:
    def __init__(self, name):
        self.start_time = datetime.now()
        self.name = name

    def elapsed(self, next_name):
        now = datetime.now()
        seconds = (now - self.start_time).total_seconds()
        report = f"{self.name} phase took {seconds:d} seconds"
        self.start_time = now
        self.name = next_name
        return report


class Sanitizer(object):
    """clean up and count potential problems with sequence

       potential problems are:
          dashes:    (optional, removed if remove_dashes=True)
          alphabet:  if not in IUPAC set, changed to 'X'
    """

    def __init__(self, remove_dashes=False):
        self.remove_dashes = remove_dashes
        self.seqs_sanitized = 0
        self.chars_in = 0
        self.chars_removed = 0
        self.chars_fixed = 0
        self.endchars_removed = 0

    def char_remover(self, s, character):
        """remove positions with a given character

        :param s: mutable sequence
        :return: sequence with characters removed
        """
        removals = [i for i, j in enumerate(s) if j == character]
        self.chars_removed += len(removals)
        [s.pop(pos - k) for k, pos in enumerate(removals)]
        return s

    def fix_alphabet(self, s):
        """replace everything out of alphabet with 'X'

        :param s: mutable sequence, upper-cased
        :return: fixed sequence
        """
        fix_positions = [pos for pos, char in enumerate(s) if char not in ALPHABET]
        self.chars_fixed = len(fix_positions)
        [s.__setitem__(pos, "X") for pos in fix_positions]
        return s

    def remove_char_on_ends(self, s, character):
        """remove leading/trailing ambiguous residues

        :param s: mutable sequence
        :return: sequence with characterss removed from ends
        """
        in_len = len(s)
        while s[-1] == character:
            s.pop()
        while s[0] == character:
            s.pop(0)
        self.endchars_removed += in_len - len(s)
        return s

    def sanitize(self, s):
        """sanitize alphabet use while checking lengths

        :param s: mutable sequence
        :return: sanitized sequence
        """
        self.seqs_sanitized += 1
        self.chars_in += len(s)
        if len(s) and self.remove_dashes:
            s = self.char_remover(s, "-")
        if len(s):
            s = self.fix_alphabet(s)
        if len(s):
            s = self.remove_char_on_ends(s, "X")
        return s


def read_synonyms(filepath):
    synonym_dict = {}
    try:
        df = pd.read_csv(filepath, sep="\t")
    except FileNotFoundError:
        logger.error(f'Synonym tsv file "{substrpath}" does not exist')
        sys.exit(1)
    except pd.errors.EmptyDataError:
        logger.error(f'Synonym tsv "{substrpath}" is empty')
        sys.exit(1)
    if len(df):
        if "#file" in df:
            df.drop("#file", axis=1, inplace=True)
        key = list(set(("Substr", "Dups")).intersection(set(df.columns)))[0]
        for group in df.groupby("id"):
            synonym_dict[group[0]] = [substr for substr in group[1][key]]
    return synonym_dict


def parse_usearch_log(filepath, rundict):
    with filepath.open() as logfile:
        for lineno, line in enumerate(logfile):
            if lineno < FIRST_LOG_LINE:
                if lineno == SEQ_IN_LINE:
                    split = line.split()
                    rundict["seqs_in"] = int(split[0])
                    rundict["singleton_seqs_in"] = int(split[4])
                if lineno == IDENT_STATS_LINE:
                    split = line.split()
                    rundict["max_identical_seqs"] = int(split[6].rstrip(","))
                    rundict["avg_identical_seqs"] = float(split[8])
                continue
            if lineno > LAST_LOG_LINE:
                break
            split = line.split()
            if split:
                stat = split[0].lower()
                if split[1] in STAT_SUFFIXES:
                    stat += "_" + split[1]
                    val = split[2]
                else:
                    val = split[1].rstrip(",")
                # rename poorly-named stats
                if stat in RENAME_STATS:
                    stat = RENAME_STATS[stat]
                # strip stats with units at the end
                conversion_factor = 1
                for unit in UNITS.keys():
                    if val.endswith(unit):
                        val = val.rstrip(unit)
                        conversion_factor = UNITS[unit]["factor"]
                        stat += "_" + UNITS[unit]["outunits"]
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


@contextlib.contextmanager
def in_working_directory(path):
    """Changes working directory and returns to previous on exit."""
    original_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original_cwd)


def get_fasta_ids(fasta):
    idset = set()
    with fasta.open() as f:
        for line in f:
            if line.startswith(">"):
                idset.add(line.split()[0][1:])
    return list(idset)


def parse_chromosome(id):
    #
    # If id contains an underscore, work on the
    # last part only (e.g., MtrunA17_Chr4g0009691)
    #
    undersplit = id.split("_")
    if len(undersplit) > 1:
        id = undersplit[-1].upper()
        if id.startswith("CHR"):
            id = id[3:]
    #
    # Chromosome numbers are integers suffixed by 'G'
    #
    try:
        chr = "Chr" + str(int(id[: id.index("G")]))
    except:
        chr = None
    return chr


def parse_subids(id):
    subids = id.split(ID_SEPARATOR)
    subids += [chr for chr in [parse_chromosome(id) for id in subids] if chr is not None]
    return subids


def parse_clusters(outdir, identity, delete=True, count_clusters=True, synonyms={}):
    cluster_list = []
    id_list = []
    degree_list = []
    size_list = []
    degree_counter = Counter()
    any_counter = Counter()
    all_counter = Counter()
    graph = nx.Graph()
    for fasta in outdir.glob("*"):
        cluster_id = int(fasta.name)
        ids = get_fasta_ids(fasta)
        if len(synonyms):
            syn_ids = set(ids).intersection(synonyms.keys())
            [ids.extend(synonyms[i]) for i in syn_ids]
        n_ids = len(ids)
        degree_list.append(n_ids)
        degree_counter.update({n_ids: 1})
        id_list += ids
        cluster_list += [cluster_id] * n_ids
        size_list += [n_ids] * n_ids
        #
        # Do 'any' and 'all' counters
        #
        id_counter = Counter()
        id_counter.update(chain.from_iterable([parse_subids(id) for id in ids]))
        if count_clusters:
            any_counter.update(id_counter.keys())
            all_counter.update([id for id in id_counter.keys() if id_counter[id] == n_ids])
        elif n_ids > 1:
            any_counter.update({s: n_ids for s in id_counter.keys()})
            all_counter.update({id: n_ids for id in id_counter.keys() if id_counter[id] == n_ids})
        #
        # Do graph components
        #
        graph.add_nodes_from(ids)
        if n_ids > 1:
            edges = combinations(ids, 2)
            graph.add_edges_from(edges, weight=n_ids)
        if delete:
            fasta.unlink()
    if delete:
        outdir.rmdir()
    return graph, cluster_list, id_list, size_list, degree_list, degree_counter, any_counter, all_counter


def prettyprint_float(x, digits):
    format_string = "%." + f"{digits:d}" + "f"
    return (format_string % x).rstrip("0").rstrip(".")


def cluster_set_name(stem, identity):
    """Get a setname that identity."""
    if identity == 1.0:
        digits = "10000"
    else:
        digits = (f"{identity:.4f}")[2:]
    return f"{stem}-nr-{digits}"


@cli.command()
@click_loguru.init_logger()
@click.argument("seqfile")
@click.option("--identity", "-i", default=0.0, help="Minimum sequence identity (float, 0-1). [default: lowest]")
@click.option("--min_id_freq", "-m", default=0, show_default=True, help="Minimum frequency of ID components.")
@click.option(
    "--delete/--no-delete",
    "-d/-n",
    is_flag=True,
    default=True,
    help="Delete primary output of clustering. [default: delete]",
)
@click.option(
    "--write_ids/--no-write_ids",
    "-w",
    is_flag=True,
    default=False,
    help="Write file of ID-to-clusters. [default: delete]",
)
@click.option(
    "--do_calc/--no-do_calc",
    "-c/-x",
    is_flag=True,
    default=True,
    help="Write file of ID-to-clusters. [default: do_calc]",
)
@click.option("--substrs", help="subpath to file of substrings. [default: none]")
@click.option("--dups", help="subpath to file of duplicates. [default: none]")
def usearch_cluster(
    seqfile, identity, delete=True, write_ids=False, do_calc=True, min_id_freq=0, substrs=None, dups=None
):
    """Cluster above a global sequence identity threshold"""
    from sh import usearch

    try:
        inpath, dirpath = get_paths_from_file(seqfile)
    except FileNotFoundError:
        logger.error(f'Input file "{seqfile}" does not exist!')
        sys.exit(1)
    stem = inpath.stem
    dirpath = inpath.parent
    outname = cluster_set_name(stem, identity)
    outdir = f"{outname}/"
    logfile = f"{outname}.log"
    outfilepath = dirpath / outdir
    logfilepath = dirpath / logfile
    histfilepath = dirpath / (f"{outname}-degreedist.tsv")
    gmlfilepath = dirpath / (f"{outname}.gml")
    statfilepath = dirpath / (f"{outname}-stats.tsv")
    anyfilepath = dirpath / (f"{outname}-anyhist.tsv")
    allfilepath = dirpath / (f"{outname}-allhist.tsv")
    idpath = dirpath / (f"{outname}-ids.tsv")
    logger.info(f"%{prettyprint_float(identity * 100, 2)}% sequence identity output to {outname}*")
    if not delete:
        logger.debug(f"Cluster files will be kept in {logfile} and {outdir}")
    if write_ids:
        logger.debug(f"File of cluster ID usage will be written to {anyfilepath} and {allfilepath}")
    if not do_calc:
        if not logfilepath.exists():
            logger.error("Previous results must exist, rerun with --do_calc")
            sys.exit(1)
        logger.debug("Using previous results for calculation")
    if min_id_freq:
        logger.debug(f"Minimum number of times ID's must occur to be counted: {min_id_freq}")
    synonyms = {}
    if substrs is not None:
        logger.debug(f"using duplicates in {dirpath / dups}")
        synonyms.update(read_synonyms(dirpath / substrs))
    if dups is not None:
        logger.debug(f"using duplicates in {dirpath/dups}")
        synonyms.update(read_synonyms(dirpath / dups))
    timer = ElapsedTimeReport("usearch")
    if do_calc:
        #
        # Delete previous results, if any.
        #
        if outfilepath.exists() and outfilepath.is_file():
            outfilepath.unlink()
        elif outfilepath.exists() and outfilepath.is_dir():
            for file in outfilepath.glob("*"):
                file.unlink()
                pass
        else:
            outfilepath.mkdir()
        #
        # Do the calculation.
        #
        with in_working_directory(dirpath):
            output = usearch(["-cluster_fast", seqfile, "-id", identity, "-clusters", outdir, "-log", logfile])
            logger.debug(output)
            logger.debug(timer.elapsed("parse"))
    run_stat_dict = OrderedDict([("divergence", 1.0 - identity)])
    parse_usearch_log(logfilepath, run_stat_dict)
    run_stats = pd.DataFrame(list(run_stat_dict.items()), columns=["stat", "val"])
    run_stats.set_index("stat", inplace=True)
    run_stats.to_csv(statfilepath, sep="\t")
    if delete:
        logfilepath.unlink()
    cluster_graph, clusters, ids, sizes, degrees, degree_counts, any_counts, all_counts = parse_clusters(
        outfilepath, identity, delete=delete, synonyms=synonyms
    )
    #
    # Write out list of clusters and ids.
    #
    id_frame = pd.DataFrame.from_dict({"id": ids, "cluster": clusters, "siz": sizes})
    id_frame.sort_values("siz", ascending=False, inplace=True)
    id_frame = id_frame.reindex(["cluster", "siz", "id",], axis=1)
    id_frame.reset_index(inplace=True)
    id_frame.drop(["index"], axis=1, inplace=True)
    id_frame.to_csv(idpath, sep="\t")
    del ids, clusters, sizes, id_frame
    logger.debug(timer.elapsed("graph"))
    #
    # Write out degree distribution.
    #
    cluster_hist = pd.DataFrame(list(degree_counts.items()), columns=["degree", "clusters"])
    # cluster_hist['degree'] = cluster_hist['degree'] - 1
    cluster_hist.sort_values(["degree"], inplace=True)
    cluster_hist.set_index("degree", inplace=True)
    total_clusters = cluster_hist["clusters"].sum()
    cluster_hist["pct_total"] = cluster_hist["clusters"] * 100.0 / total_clusters
    cluster_hist.to_csv(histfilepath, sep="\t", float_format="%06.3f")
    del degree_counts
    #
    # Do histograms of "any" and "all" id usage in cluster
    #
    hist_value = f"{identity:f}"
    any_hist = pd.DataFrame(list(any_counts.items()), columns=["id", hist_value])
    any_hist.set_index("id", inplace=True)
    any_hist.sort_values(hist_value, inplace=True, ascending=False)
    all_hist = pd.DataFrame(list(all_counts.items()), columns=["id", hist_value])
    all_hist.set_index("id", inplace=True)
    all_hist.sort_values(hist_value, inplace=True, ascending=False)
    if min_id_freq:
        any_hist = any_hist[any_hist[hist_value] > min_id_freq]
        all_hist = all_hist[all_hist[hist_value] > min_id_freq]
    if write_ids:
        any_hist.to_csv(anyfilepath, sep="\t")
        all_hist.to_csv(allfilepath, sep="\t")
    #
    # Compute cluster stats
    #
    # degree_sequence = sorted([d for n, d in cluster_graph.degree()], reverse=True)
    # degreeCount = Counter(degree_sequence)
    # degree_hist = pd.DataFrame(list(degreeCount.items()),
    #                           columns=['degree', 'count'])
    # degree_hist.set_index('degree', inplace=True)
    # degree_hist.sort_values('degree', inplace=True)
    # degree_hist.to_csv(histfilepath, sep='\t')
    nx.write_gml(cluster_graph, gmlfilepath)
    logger.debug(timer.elapsed("final"))
    return run_stat_dict, cluster_graph, cluster_hist, any_hist, all_hist


@cli.command()
@click_loguru.init_logger()
@click.argument("seqfile")
@click.option("--steps", "-s", default=DEFAULT_STEPS, show_default=True, help="# of steps from lowest to highest")
@click.option("--min_id_freq", "-m", default=0, show_default=True, help="Minimum frequency of ID components.")
@click.option("--substrs", help="subpath to file of substrings. [default: none]")
@click.option("--dups", help="subpath to file of duplicates. [default: none]")
def cluster_in_steps(seqfile, steps, min_id_freq=0, substrs=None, dups=None):
    """Cluster in steps from low to 100% sequence identity"""
    try:
        inpath, dirpath = get_paths_from_file(seqfile)
    except FileNotFoundError:
        logger.error('Input file "%s" does not exist!', seqfile)
        sys.exit(1)
    stat_path = dirpath / (inpath.stem + STATFILE_SUFFIX)
    any_path = dirpath / (inpath.stem + ANYFILE_SUFFIX)
    all_path = dirpath / (inpath.stem + ALLFILE_SUFFIX)
    logsteps = [1.0] + list(1.0 - np.logspace(IDENT_LOG_MIN, IDENT_LOG_MAX, num=steps))
    logger.info(
        "Clustering at %d levels from %s%% to %s%% global sequence identity",
        steps,
        prettyprint_float(min(logsteps) * 100.0, 2),
        prettyprint_float(max(logsteps) * 100.0, 2),
    )
    stat_list = []
    all_frames = []
    any_frames = []
    for id_level in logsteps:
        stats, graph, hist, any, all = usearch_cluster.callback(
            seqfile, id_level, min_id_freq=min_id_freq, substrs=substrs, dups=dups
        )
        stat_list.append(stats)
        any_frames.append(any)
        all_frames.append(all)
    logger.info("Collating results on %s", seqfile)
    #
    # Concatenate and write stats
    #
    stats = pd.DataFrame(stat_list)
    stats.to_csv(stat_path, sep="\t")
    #
    # Concatenate any/all data
    #
    any = pd.concat(any_frames, axis=1, join="inner", sort=True, ignore_index=False)
    any.to_csv(any_path, sep="\t")
    all = pd.concat(all_frames, axis=1, join="inner", sort=True, ignore_index=False)
    all.to_csv(all_path, sep="\t")


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("infile")
def clusters_to_histograms(infile):
    """Compute histograms from a tab-delimited cluster file"""
    try:
        inpath, dirpath = get_paths_from_file(infile)
    except FileNotFoundError:
        logger.error('Input file "%s" does not exist!', infile)
        sys.exit(1)
    histfilepath = dirpath / (inpath.stem + "-sizedist.tsv")
    clusters = pd.read_csv(dirpath / infile, sep="\t", index_col=0)
    cluster_counter = Counter()
    for cluster_id, group in clusters.groupby(["cluster"]):
        cluster_counter.update({len(group): 1})
    logger.info("writing to %s", histfilepath)
    cluster_hist = pd.DataFrame(list(cluster_counter.items()), columns=["siz", "clusts"])
    total_clusters = cluster_hist["clusts"].sum()
    cluster_hist["%clusts"] = cluster_hist["clusts"] * 100.0 / total_clusters
    cluster_hist["%genes"] = cluster_hist["clusts"] * cluster_hist["siz"] * 100.0 / len(clusters)
    cluster_hist.sort_values(["siz"], inplace=True)
    cluster_hist.set_index("siz", inplace=True)
    cluster_hist.to_csv(histfilepath, sep="\t", float_format="%06.3f")


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("file1")
@click.argument("file2")
def compare_clusters(file1, file2):
    """ compare one cluster file with another
    """
    path1 = Path(file1)
    path2 = Path(file2)
    commondir = Path(os.path.commonpath([path1, path2]))
    missing1 = commondir / "notin1.tsv"
    missing2 = commondir / "notin2.tsv"
    clusters1 = pd.read_csv(path1, sep="\t", index_col=0)
    print(f"{len(clusters1):d} members in {file1}")
    clusters2 = pd.read_csv(path2, sep="\t", index_col=0)
    print(f"{len(clusters2):d} members in {file2}")
    ids1 = set(clusters1["id"])
    ids2 = set(clusters2["id"])
    notin1 = pd.DataFrame(ids2.difference(ids1), columns=["id"])
    notin1.sort_values("id", inplace=True)
    notin1.to_csv(missing1, sep="\t")
    notin2 = pd.DataFrame(ids1.difference(ids2), columns=["id"])
    notin2.sort_values("id", inplace=True)
    notin2.to_csv(missing2, sep="\t")

    print("%d ids not in ids1" % len(notin1))
    print("%d ids not in ids2" % len(notin2))
    print(f"{len(clusters1):d} in {file1} after dropping")


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("setname")
@click.argument("filelist", nargs=-1)
def scanfiles(setname, filelist):
    """scan a set of files, producing summary statistics"""
    setpath = Path(setname)
    if len(filelist) < 1:
        logger.error("Empty FILELIST, aborting")
        sys.exit(1)
    if setpath.exists() and setpath.is_file():
        setpath.unlink()
    elif setpath.exists() and setpath.is_dir():
        logger.debug("set path %s exists", setname)
    else:
        logger.info('creating output directory "%s/"', setname)
        setpath.mkdir()
    outfile = setpath / "all.faa"
    statfile = setpath / "records.tsv"
    out_sequences = []
    lengths = []
    ids = []
    files = []
    positions = []
    for file in filelist:
        position = 0
        logger.info("scanning %s", file)
        sanitizer = Sanitizer(remove_dashes=True)
        with open(file, "rU") as handle:
            for record in SeqIO.parse(handle, SEQ_FILE_TYPE):
                seq = record.seq.upper().tomutable()
                seq = sanitizer.sanitize(seq)
                if not len(seq):
                    # zero-length string after sanitizing
                    continue
                record.seq = seq.toseq()
                ids.append(record.id)
                lengths.append(len(record))
                out_sequences.append(record)
                files.append(file)
                positions.append(position)
                position += 1
    logger.debug("writing output files")
    with outfile.open("w") as output_handle:
        SeqIO.write(out_sequences, output_handle, SEQ_FILE_TYPE)
    df = pd.DataFrame(list(zip(files, ids, positions, lengths)), columns=["file", "id", "pos", "len"])
    df.to_csv(statfile, sep="\t")


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("infile")
@click.argument("recordfile")
def add_singletons(infile, recordfile):
    """Add singleton clusters to cluster file"""
    clusters = pd.read_csv(infile, header=None, names=["cluster", "id"], sep="\t")
    records = pd.read_csv(recordfile, sep="\t")
    id_set = set(records["id"])
    records.drop(["file", "pos"], axis=1, inplace=True)
    records.set_index("id", inplace=True)
    records.drop(["Unnamed: 0"], axis=1, inplace=True)
    len_dict = records.to_dict("id")
    del records
    outfile = "clusters.tsv"
    ids = []
    cluster_ids = []
    sizes = []
    lengths = []
    n_clusters = 0
    grouping = clusters.groupby("cluster").size().sort_values(ascending=False)
    for idx in grouping.index:
        size = grouping.loc[idx]
        cluster = clusters.loc[(clusters["cluster"] == idx)]
        for gene_id in cluster["id"]:
            id_set.remove(gene_id)
            ids.append(gene_id)
            cluster_ids.append(n_clusters)
            sizes.append(size)
            lengths.append(len_dict[gene_id]["len"])
        n_clusters += 1
    logger.info("%s singletons", len(id_set))
    for singleton in id_set:
        cluster_ids.append(n_clusters)
        ids.append(singleton)
        sizes.append(1)
        lengths.append(len_dict[singleton]["len"])
        n_clusters += 1
    df = pd.DataFrame(list(zip(cluster_ids, ids, sizes, lengths)), columns=["cluster", "id", "siz", "len"])
    df.to_csv(outfile, sep="\t")


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.argument("infile")
def adjacency_to_graph(infile):
    gmlfilepath = "synteny.gml"
    clusters = pd.read_csv(infile, sep="\t", index_col=0)
    graph = nx.Graph()
    grouping = clusters.groupby("cluster").size().sort_values(ascending=False)
    for idx in grouping.index:
        size = grouping.loc[idx]
        cluster = clusters.loc[(clusters["cluster"] == idx)]
        ids = list(cluster["id"])
        #
        # Do graph components
        #
        graph.add_nodes_from(ids)
        if len(ids) > 1:
            edges = combinations(ids, 2)
            graph.add_edges_from(edges, weight=len(ids))
    nx.write_gml(graph, gmlfilepath)


def compute_subclusters(cluster, cluster_size_dict=None):
    #
    # compute a dictionary of per-subcluster stats first,
    # so subclusters can be ordered by length
    #
    subcl_frame = pd.DataFrame(
        [
            {"homo_id": i, "mean_len": g["len"].mean(), "std": g["len"].std(), "sub_siz": len(g)}
            for i, g in cluster.groupby("link")
        ]
    )
    subcl_frame["link"] = subcl_frame["homo_id"]
    subcl_frame["cont"] = subcl_frame["sub_siz"] == [cluster_size_dict[id] for id in subcl_frame["homo_id"]]
    subcl_frame.loc[subcl_frame["cont"], "link"] = np.nan
    subcl_frame.sort_values(["sub_siz", "std", "mean_len"], ascending=[False, True, False], inplace=True)
    subcl_frame["sub"] = list(range(len(subcl_frame)))
    subcl_frame.index = list(range(len(subcl_frame)))
    # normalized length is NaN for first element
    subcl_frame["norm"] = [np.nan] + list(subcl_frame["mean_len"][1:] / subcl_frame["mean_len"][0])
    subcl_dict = subcl_frame.set_index("homo_id").to_dict("index")
    #
    # subcluster attributes are done, now copy them into the cluster frame
    #
    for attr in ["norm", "std", "sub_siz", "sub", "link"]:
        cluster[attr] = [subcl_dict[id][attr] for id in cluster["link"]]
    if subcl_frame["cont"].all():
        cluster["cont"] = 1
    else:
        cluster["cont"] = 0
    return cluster


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option("--first_n", default=0, show_default=True, help="Only process this many clusters.")
@click.option("--clust_size", default=0, show_default=True, help="Process only clusters of this size.")
@click.option("--parallel/--no-parallel", is_flag=True, default=True, show_default=True, help="Process in parallel.")
@click.option("-q", "--quiet", is_flag=True, show_default=True, default=False, help="Suppress logging to stderr.")
@click.argument("synfile")
@click.argument("homofile")
def combine_clusters(first_n, clust_size, synfile, homofile, quiet, parallel):
    """combine synteny and homology clusters"""
    timer = ElapsedTimeReport("reading/preparing")
    syn = pd.read_csv(synfile, sep="\t", index_col=0)
    homo = pd.read_csv(homofile, sep="\t", index_col=0)
    homo_id_dict = dict(zip(homo.id, homo.cluster))
    cluster_size_dict = dict(zip(homo.cluster, homo.siz))
    cluster_count = 0
    arg_list = []
    if first_n:
        logger.debug(f"processing only first {first_n:d} clusters")
    if not quiet and clust_size:
        logger.info("only clusters of size %d will be used", clust_size)
    syn["link"] = [homo_id_dict[id] for id in syn["id"]]
    for cluster_id, gr in syn.groupby(["cluster"]):
        cl = gr.copy()  # copy, so mutable
        clsize = cl["siz"].iloc[0]
        if clust_size and (clsize != clust_size):
            continue
        if first_n and cluster_count == first_n:
            break
        arg_list.append(cl)
        cluster_count += 1
    del syn
    if parallel:
        bag = db.from_sequence(arg_list)
    else:
        cluster_list = []
    if not quiet:
        logger.info("Combining %d synteny/homology clusters:", cluster_count)
        ProgressBar().register()
    if parallel:
        cluster_list = bag.map(compute_subclusters, cluster_size_dict=cluster_size_dict)
    else:
        for clust in arg_list:
            cluster_list.append(compute_subclusters(clust, cluster_size_dict))
    logger.debug(timer.elapsed("concatenating/writing clusters"))
    out_frame = pd.concat(cluster_list)
    out_frame.sort_values(["cluster", "sub"], inplace=True)
    out_frame.index = list(range(len(out_frame)))
    for column in ["sub", "sub_siz"]:
        out_frame[column] = out_frame[column].astype(int)
    out_frame["link"] = out_frame["link"].map(lambda x: "" if pd.isnull(x) else f"{x:.0f}")
    out_frame.to_csv(
        "combined.tsv",
        sep="\t",
        columns=["cluster", "siz", "sub", "sub_siz", "cont", "norm", "std", "len", "link", "id"],
        float_format="%.3f",
        index=False,
    )
    logger.debug(timer.elapsed("computing stats"))
    n_fully_contained = len(set(out_frame[out_frame["cont"] == 1]["cluster"]))
    logger.info(
        "%d of %d clusters are fully contained (%.1f%%)",
        n_fully_contained,
        cluster_count,
        n_fully_contained * 100 / cluster_count,
    )
    # print('subclusters\tN\t%clusts')
    # for count in sorted(contained_subclst_count.keys()):
    #    print('%d\t%d\t%.1f'%(count, contained_subclst_count[count],
    #                         contained_subclst_count[count]*100/cluster_count))
    # uncontained = cluster_count - n_fully_contained
    # print('%d of %d clusters are unontained (%.1f%%)'
    #      % (uncontained, cluster_count,
    #         uncontained * 100 / cluster_count))
