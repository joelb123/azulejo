# -*- coding: utf-8 -*-
"""Core logic for computing subtrees."""

# standard library imports
import contextlib
import os
import sys
from collections import Counter
from collections import OrderedDict
from itertools import chain
from itertools import combinations
from pathlib import Path

# third-party imports
import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO

# first-party imports
import sh

# module imports
from .common import NAME
from .common import cluster_set_name
from .common import fasta_records
from .common import get_paths_from_file
from .common import homo_degree_dist_filename
from .common import logger
from .common import protein_properties_filename
from .common import write_tsv_or_parquet
from .protein import Sanitizer

# global constants
STATFILE_SUFFIX = f"-{NAME}_stats.tsv"
ANYFILE_SUFFIX = f"-{NAME}_ids-any.tsv"
ALLFILE_SUFFIX = f"-{NAME}_ids-all.tsv"
CLUSTFILE_SUFFIX = f"-{NAME}_clusts.tsv"
SEQ_FILE_TYPE = "fasta"
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
STAT_SUFFIXES = ["size", "mem", "time", "memory"]
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
FASTA_EXT_LIST = [".faa", ".fa", ".fasta"]
FAA_EXT = "faa"


# helper functions
def read_synonyms(filepath):
    """Read a file of synonymous IDs into a dictionary."""
    synonym_dict = {}
    try:
        synonym_frame = pd.read_csv(filepath, sep="\t")
    except FileNotFoundError:
        logger.error(f'Synonym tsv file "{filepath}" does not exist')
        sys.exit(1)
    except pd.errors.EmptyDataError:
        logger.error(f'Synonym tsv "{filepath}" is empty')
        sys.exit(1)
    if len(synonym_frame) > 0:
        if "#file" in synonym_frame:
            synonym_frame.drop("#file", axis=1, inplace=True)
        key = list({("Substr", "Dups")}.intersection({synonym_frame.columns}))[
            0
        ]
        for group in synonym_frame.groupby("id"):
            synonym_dict[group[0]] = group[1][key]
    return synonym_dict


def parse_usearch_log(filepath, rundict):
    """Parse the usearch log file into a stats dictionary."""
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
                stat = RENAME_STATS.get(stat, stat)
                # strip stats with units at the end
                conversion_factor = 1
                for unit in UNITS:
                    if val.endswith(unit):
                        val = val.rstrip(unit)
                        conversion_factor = UNITS[unit]["factor"]
                        stat += "_" + UNITS[unit]["outunits"]
                        break
                # convert string values to int or float where possible
                try:
                    val = int(val)
                    val *= conversion_factor
                except ValueError:
                    try:
                        val = float(val)
                        val *= conversion_factor
                    except ValueError:
                        pass
                rundict[stat] = val


@contextlib.contextmanager
def in_working_directory(path):
    """Change working directory and return to previous wd on exit."""
    original_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original_cwd)


def get_fasta_ids(fasta):
    """Get the IDS from a FASTA file."""
    idset = set()
    with fasta.open() as fasta_fh:
        for line in fasta_fh:
            if line.startswith(">"):
                idset.add(line.split()[0][1:])
    return list(idset)


def parse_chromosome(ident):
    """Parse chromosome identifiers."""
    # If ident contains an underscore, work on the
    # last part only (e.g., MtrunA17_Chr4g0009691)
    undersplit = ident.split("_")
    if len(undersplit) > 1:
        ident = undersplit[-1].upper()
        if ident.startswith("CHR"):
            ident = ident[3:]
    # Chromosome numbers are integers suffixed by 'G'
    try:
        chromosome = "Chr" + str(int(ident[: ident.index("G")]))
    except ValueError:
        chromosome = None
    return chromosome


def parse_subids(ident):
    """Parse the subidentifiers from identifiers."""
    subids = ident.split(ID_SEPARATOR)
    subids += [
        chromosome
        for chromosome in [parse_chromosome(ident) for ident in subids]
        if chromosome is not None
    ]
    return subids


def parse_clusters(outdir, delete=True, count_clusters=True, synonyms=None):
    """Parse clusters, counting occurrances."""
    if synonyms is None:
        synonyms = {}
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
        if len(synonyms) > 0:
            syn_ids = set(ids).intersection(synonyms.keys())
            for i in syn_ids:
                ids.extend(synonyms[i])
        n_ids = len(ids)
        degree_list.append(n_ids)
        degree_counter.update({n_ids: 1})
        id_list += ids
        cluster_list += [cluster_id] * n_ids
        size_list += [n_ids] * n_ids
        # Do 'any' and 'all' counters
        id_counter = Counter()
        id_counter.update(
            chain.from_iterable(
                [parse_subids(cluster_id) for cluster_id in ids]
            )
        )
        if count_clusters:
            any_counter.update(id_counter.keys())
            all_counter.update(
                [
                    cluster_id
                    for cluster_id in id_counter.keys()
                    if id_counter[cluster_id] == n_ids
                ]
            )
        elif n_ids > 1:
            any_counter.update({s: n_ids for s in id_counter.keys()})
            all_counter.update(
                {
                    cluster_id: n_ids
                    for cluster_id in id_counter.keys()
                    if id_counter[cluster_id] == n_ids
                }
            )
        # Do graph components
        graph.add_nodes_from(ids)
        if n_ids > 1:
            edges = combinations(ids, 2)
            graph.add_edges_from(edges, weight=n_ids)
        if delete:
            fasta.unlink()
    if delete:
        outdir.rmdir()
    return (
        graph,
        cluster_list,
        id_list,
        size_list,
        degree_list,
        degree_counter,
        any_counter,
        all_counter,
    )


def prettyprint_float(val, digits):
    """Print a floating-point value in a nice way."""
    format_string = "%." + f"{digits:d}" + "f"
    return (format_string % val).rstrip("0").rstrip(".")


def homology_cluster(
    seqfile,
    identity,
    delete=True,
    write_ids=False,
    do_calc=True,
    min_id_freq=0,
    substrs=None,
    dups=None,
    cluster_stats=True,
    outname=None,
    click_loguru=None,
):
    """Cluster at a global sequence identity threshold."""
    try:
        usearch = sh.Command("usearch")
    except sh.CommandNotFound:
        logger.error("usearch must be installed first.")
        sys.exit(1)
    try:
        inpath, dirpath = get_paths_from_file(seqfile)
    except FileNotFoundError:
        logger.error(f'Input file "{seqfile}" does not exist!')
        sys.exit(1)
    stem = inpath.stem
    dirpath = inpath.parent
    if outname is None:
        outname = cluster_set_name(stem, identity)
    outdir = f"{outname}/"
    logfile = f"{outname}.log"
    outfilepath = dirpath / outdir
    logfilepath = dirpath / logfile
    histfilepath = dirpath / homo_degree_dist_filename(outname)
    gmlfilepath = dirpath / f"{outname}.gml"
    statfilepath = dirpath / f"{outname}-stats.tsv"
    anyfilepath = dirpath / f"{outname}-anyhist.tsv"
    allfilepath = dirpath / f"{outname}-allhist.tsv"
    idpath = dirpath / f"{outname}-ids.tsv"
    if identity == 0.0:
        identity_string = "Minimum"
    else:
        identity_string = f"{prettyprint_float(identity *100, 2)}%"
    logger.info(f'{identity_string} sequence identity cluster "{outname}":')
    if not delete:
        logger.debug(f"Cluster files will be kept in {logfile} and {outdir}")
    if cluster_stats and write_ids:
        logger.debug(
            f"File of cluster ID usage will be written to {anyfilepath} and"
            f" {allfilepath}"
        )
    if not do_calc:
        if not logfilepath.exists():
            logger.error("Previous results must exist, rerun with --do_calc")
            sys.exit(1)
        logger.debug("Using previous results for calculation")
    if min_id_freq:
        logger.debug(
            "Minimum number of times ID's must occur to be counted:"
            f" {min_id_freq}"
        )
    synonyms = {}
    if substrs is not None:
        logger.debug(f"using duplicates in {dirpath / dups}")
        synonyms.update(read_synonyms(dirpath / substrs))
    if dups is not None:
        logger.debug(f"using duplicates in {dirpath/dups}")
        synonyms.update(read_synonyms(dirpath / dups))
    click_loguru.elapsed_time("Clustering")
    if do_calc:
        #
        # Delete previous results, if any.
        #
        if outfilepath.exists() and outfilepath.is_file():
            outfilepath.unlink()
        elif outfilepath.exists() and outfilepath.is_dir():
            for file in outfilepath.glob("*"):
                file.unlink()
        else:
            outfilepath.mkdir()
        #
        # Do the calculation.
        #
        with in_working_directory(dirpath):
            output = usearch(
                [
                    "-cluster_fast",
                    seqfile,
                    "-id",
                    identity,
                    "-clusters",
                    outdir,
                    "-log",
                    logfile,
                ]
            )
            logger.debug(output)
    run_stat_dict = OrderedDict([("divergence", 1.0 - identity)])
    parse_usearch_log(logfilepath, run_stat_dict)
    run_stats = pd.DataFrame(
        list(run_stat_dict.items()), columns=["stat", "val"]
    )
    run_stats.set_index("stat", inplace=True)
    write_tsv_or_parquet(run_stats, statfilepath)
    if delete:
        logfilepath.unlink()
    if not cluster_stats:
        file_sizes = []
        file_names = []
        record_counts = []
        logger.debug("Ordering clusters by number of records and size.")
        for fasta_path in outfilepath.glob("*"):
            records, size = fasta_records(fasta_path)
            if records == 1:
                fasta_path.unlink()
                continue
            file_names.append(fasta_path.name)
            file_sizes.append(size)
            record_counts.append(records)
        file_frame = pd.DataFrame(
            list(zip(file_names, file_sizes, record_counts)),
            columns=["name", "size", "seqs"],
        )
        file_frame.sort_values(
            by=["seqs", "size"], ascending=False, inplace=True
        )
        file_frame["idx"] = range(len(file_frame))
        for unused_id, row in file_frame.iterrows():
            (outfilepath / row["name"]).rename(
                outfilepath / f'{row["idx"]}.fa'
            )
        file_frame.drop(["name"], axis=1, inplace=True)
        file_frame.set_index("idx", inplace=True)
        # write_tsv_or_parquet(file_frame, "clusters.tsv")
        # cluster histogram
        cluster_hist = pd.DataFrame(file_frame["seqs"].value_counts())
        cluster_hist.rename(columns={"seqs": "clusters"}, inplace=True)
        cluster_hist.index.name = "n"
        cluster_hist.sort_index(inplace=True)
        total_seqs = sum(file_frame["seqs"])
        n_clusters = len(file_frame)
        cluster_hist["pct_clusts"] = (
            cluster_hist["clusters"] * 100.0 / n_clusters
        )
        cluster_hist["pct_seqs"] = (
            cluster_hist["clusters"] * cluster_hist.index * 100.0 / total_seqs
        )
        cluster_hist.to_csv(
            "cluster_hist.tsv", sep="\t", float_format="%06.3f"
        )
        return n_clusters, run_stats, cluster_hist
    (
        cluster_graph,
        clusters,
        ids,
        sizes,
        unused_degrees,
        degree_counts,
        any_counts,
        all_counts,
    ) = parse_clusters(  # pylint: disable=unused-variable
        outfilepath, delete=delete, synonyms=synonyms
    )
    #
    # Write out list of clusters and ids.
    #
    id_frame = pd.DataFrame.from_dict(
        {
            "id": ids,
            "hom.cluster": pd.array(clusters, dtype=pd.UInt32Dtype()),
            "siz": sizes,
        }
    )
    id_frame.sort_values("siz", ascending=False, inplace=True)
    id_frame = id_frame.reindex(["hom.cluster", "siz", "id"], axis=1,)
    id_frame.reset_index(inplace=True)
    id_frame.drop(["index"], axis=1, inplace=True)
    id_frame.to_csv(idpath, sep="\t")
    del ids, clusters, sizes, id_frame
    click_loguru.elapsed_time("graph")
    #
    # Write out degree distribution.
    #
    cluster_hist = pd.DataFrame(
        list(degree_counts.items()), columns=["degree", "clusters"]
    )
    cluster_hist.sort_values(["degree"], inplace=True)
    cluster_hist.set_index("degree", inplace=True)
    total_clusters = cluster_hist["clusters"].sum()
    cluster_hist["pct_total"] = (
        cluster_hist["clusters"] * 100.0 / total_clusters
    )
    cluster_hist.to_csv(histfilepath, sep="\t", float_format="%06.3f")
    del degree_counts
    #
    # Do histograms of "any" and "all" id usage in cluster
    #
    hist_value = f"{identity:f}"
    any_hist = pd.DataFrame(
        list(any_counts.items()), columns=["id", hist_value]
    )
    any_hist.set_index("id", inplace=True)
    any_hist.sort_values(hist_value, inplace=True, ascending=False)
    all_hist = pd.DataFrame(
        list(all_counts.items()), columns=["id", hist_value]
    )
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
    click_loguru.elapsed_time("final")
    return run_stats, cluster_graph, cluster_hist, any_hist, all_hist


def cluster_in_steps(seqfile, steps, min_id_freq=0, substrs=None, dups=None):
    """Cluster in steps from low to 100% identity."""
    try:
        inpath, dirpath = get_paths_from_file(seqfile)
    except FileNotFoundError:
        logger.error('Input file "%s" does not exist!', seqfile)
        sys.exit(1)
    stat_path = dirpath / (inpath.stem + STATFILE_SUFFIX)
    any_path = dirpath / (inpath.stem + ANYFILE_SUFFIX)
    all_path = dirpath / (inpath.stem + ALLFILE_SUFFIX)
    logsteps = [1.0] + list(
        1.0 - np.logspace(IDENT_LOG_MIN, IDENT_LOG_MAX, num=steps)
    )
    min_fmt = prettyprint_float(min(logsteps) * 100.0, 2)
    max_fmt = prettyprint_float(max(logsteps) * 100.0, 2)
    logger.info(
        f"Clustering at {steps} levels from {min_fmt}% to {max_fmt}% global"
        " sequence identity"
    )
    stat_list = []
    all_frames = []
    any_frames = []
    for id_level in logsteps:
        (
            stats,
            unused_graph,
            unused_hist,
            any_,
            all_,
        ) = homology_cluster(  # pylint: disable=unused-variable
            seqfile,
            id_level,
            min_id_freq=min_id_freq,
            substrs=substrs,
            dups=dups,
        )
        stat_list.append(stats)
        any_frames.append(any_)
        all_frames.append(all_)
    logger.info(f"Collating results on {seqfile}.")
    #
    # Concatenate and write stats
    #
    stats = pd.DataFrame(stat_list)
    stats.to_csv(stat_path, sep="\t")
    #
    # Concatenate any/all data
    #
    any_ = pd.concat(
        any_frames, axis=1, join="inner", sort=True, ignore_index=False
    )
    any_.to_csv(any_path, sep="\t")
    all_ = pd.concat(
        all_frames, axis=1, join="inner", sort=True, ignore_index=False
    )
    all_.to_csv(all_path, sep="\t")


def clusters_to_histograms(infile):
    """Compute histograms from cluster file."""
    try:
        inpath, dirpath = get_paths_from_file(infile)
    except FileNotFoundError:
        logger.error(f'Input file "{infile}" does not exist!')
        sys.exit(1)
    histfilepath = dirpath / f"{inpath.stem}-sizedist.tsv"
    clusters = pd.read_csv(dirpath / infile, sep="\t", index_col=0)
    cluster_counter = Counter()
    for unused_cluster_id, group in clusters.groupby(
        ["hom.cluster"]
    ):  # pylint: disable=unused-variable
        cluster_counter.update({len(group): 1})
    logger.info(f"Writing to {histfilepath}.")
    cluster_hist = pd.DataFrame(
        list(cluster_counter.items()), columns=["siz", "clusts"]
    )
    total_clusters = cluster_hist["clusts"].sum()
    cluster_hist["%clusts"] = cluster_hist["clusts"] * 100.0 / total_clusters
    cluster_hist["%genes"] = (
        cluster_hist["clusts"] * cluster_hist["siz"] * 100.0 / len(clusters)
    )
    cluster_hist.sort_values(["siz"], inplace=True)
    cluster_hist.set_index("siz", inplace=True)
    cluster_hist.to_csv(histfilepath, sep="\t", float_format="%06.3f")


def compare_clusters(file1, file2):
    """Compare cluster files."""
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


def cleanup_fasta(
    set_path,
    fasta_path_or_handle,
    filestem,
    write_fasta=True,
    write_stats=True,
):
    """Sanitize and characterize protein FASTA files."""
    out_sequences = []
    ids = []
    lengths = []
    n_ambigs = []
    m_starts = []
    no_stops = []
    seqs = []
    sanitizer = Sanitizer(filestem)
    if hasattr(fasta_path_or_handle, "read"):
        handle = fasta_path_or_handle
    else:
        handle = fasta_path_or_handle.open("rU")
    with handle:
        for record in SeqIO.parse(handle, SEQ_FILE_TYPE):
            seq = record.seq.upper().tomutable()
            try:
                seq, m_start, no_stop, n_ambig = sanitizer.sanitize(seq)
            except ValueError:  # zero-length sequence after sanitizing
                continue
            seqs.append(str(seq))
            if write_fasta:
                record.seq = seq.toseq()
                out_sequences.append(record)
            ids.append(record.id)
            lengths.append(len(seq))
            n_ambigs.append(n_ambig)
            m_starts.append(m_start)
            no_stops.append(no_stop)
    if write_fasta:
        with (set_path / f"{filestem}.fa").open("w") as output_handle:
            SeqIO.write(out_sequences, output_handle, SEQ_FILE_TYPE)
    properties_frame = pd.DataFrame(
        list(zip(ids, lengths, n_ambigs, m_starts, no_stops, seqs)),
        columns=[
            "ID",
            "prot.len",
            "prot.n_ambig",
            "prot.m_start",
            "prot.no_stop",
            "prot.seq",
        ],
    )
    properties_frame = properties_frame.set_index(["ID"])
    if write_stats:
        propfile_path = set_path / protein_properties_filename(filestem)
        properties_frame.to_csv(propfile_path, sep="\t")
        frame_return = len(properties_frame)
    else:
        propfile_path = None
        frame_return = properties_frame
    return (
        filestem,
        propfile_path,
        frame_return,
        sanitizer.file_stats(),
    )


def compute_subclusters(cluster, cluster_size_dict=None):
    """Compute dictionary of per-subcluster stats."""
    subcl_frame = pd.DataFrame(
        [
            {
                "homo_id": i,
                "mean_len": g["len"].mean(),
                "std": g["len"].std(),
                "sub_siz": len(g),
            }
            for i, g in cluster.groupby("link")
        ]
    )
    subcl_frame["link"] = subcl_frame["homo_id"]
    subcl_frame["cont"] = subcl_frame["sub_siz"] == [
        cluster_size_dict[cluster_id] for cluster_id in subcl_frame["homo_id"]
    ]
    subcl_frame.loc[subcl_frame["cont"], "link"] = pd.NA
    subcl_frame.sort_values(
        ["sub_siz", "std", "mean_len"],
        ascending=[False, True, False],
        inplace=True,
    )
    subcl_frame["sub"] = list(range(len(subcl_frame)))
    subcl_frame.index = list(range(len(subcl_frame)))
    # normalized length is NaN for first element
    subcl_frame["norm"] = [pd.NA] + list(
        subcl_frame["mean_len"][1:] / subcl_frame["mean_len"][0]
    )
    subcl_dict = subcl_frame.set_index("homo_id").to_dict("index")
    #
    # subcluster attributes are done, now copy them into the cluster frame
    #
    for attr in ["norm", "std", "sub_siz", "sub", "link"]:
        cluster[attr] = [
            subcl_dict[cluster_id][attr] for cluster_id in cluster["link"]
        ]
    if subcl_frame["cont"].all():
        cluster["cont"] = 1
    else:
        cluster["cont"] = 0
    return cluster
