# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import os
import statistics
import sys
from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import click
import gffpandas.gffpandas as gffpd
import numpy as np
import pandas as pd
import sh
from Bio import SeqIO
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import FAA_EXT
from .common import GFF_EXT
from .core import cluster_set_name
from .core import prepare_protein_files
from .core import usearch_cluster

# global constants


HOMOLOGY_ENDING = "-homology.tsv"
FILES_ENDING = "-files.tsv"
SYNTENY_ENDING = "-synteny.tsv"
PROXY_ENDING = "-proxy.tsv"


def synteny_block_func(k, rmer, frame, name_only=False):
    """Return a synteny block closure and its name."""
    if name_only and rmer:
        return f"rmer{k}"
    if name_only and not rmer:
        return f"kmer{k}"
    frame_len = len(frame)
    cluster_size_col = frame.columns.get_loc("cluster_size")
    cluster_col = frame.columns.get_loc("cluster_id")

    def kmer_block(first_index):
        """Calculate a reversible hash of cluster values.."""
        cluster_list = []
        for idx in range(first_index, first_index + k):
            if idx + 1 > frame_len or frame.iloc[idx, cluster_size_col] == 1:
                return (
                    0,
                    0,
                    0,
                )
            cluster_list.append(frame.iloc[idx, cluster_col])
        fwd_hash = hash(tuple(cluster_list))
        rev_hash = hash(tuple(reversed(cluster_list)))
        if fwd_hash > rev_hash:
            return k, 1, fwd_hash
        return k, -1, rev_hash

    def rmer_block(first_index):
        """Calculate a reversible cluster hash, ignoring repeats."""
        cluster_list = []
        idx = first_index
        last_cluster = None
        while len(cluster_list) < k:
            if idx + 1 > frame_len or frame.iloc[idx, cluster_size_col] == 1:
                return (
                    0,
                    0,
                    0,
                )
            current_cluster = frame.iloc[idx, cluster_col]
            if current_cluster == last_cluster:
                idx += 1
            else:
                last_cluster = current_cluster
                cluster_list.append(current_cluster)
                idx += 1
        fwd_hash = hash(tuple(cluster_list))
        rev_hash = hash(tuple(reversed(cluster_list)))
        if fwd_hash > rev_hash:
            return idx - first_index, 1, fwd_hash
        return idx - first_index, -1, rev_hash

    if rmer:
        return rmer_block
    return kmer_block


def read_files(setname, synteny=None):
    """Read previously-calculated homology/synteny files and file frame."""
    set_path = Path(setname)
    files_frame_path = set_path / f"{setname}{FILES_ENDING}"
    try:
        file_frame = pd.read_csv(files_frame_path, index_col=0, sep="\t")
    except FileNotFoundError:
        logger.error(f"Unable to read files frome from {files_frame_path}")
        sys.exit(1)
    if synteny is None:
        ending = HOMOLOGY_ENDING
        file_type = "homology"
    else:
        ending = f"-{synteny}{SYNTENY_ENDING}"
        file_type = "synteny"
    paths = list(set_path.glob("*" + ending))
    stems = [p.name[: -len(ending)] for p in paths]
    if len(stems) != len(file_frame):
        logger.error(
            f"Number of {file_type} files ({len(stems)})"
            + f"is not the same as length of file frame({len(file_frame)})."
        )
        sys.exit(1)
    frame_dict = {}
    for i, path in enumerate(paths):
        logger.debug(f"Reading homology file {path}.")
        frame_dict[stems[i]] = pd.read_csv(path, index_col=0, sep="\t")
    return file_frame, frame_dict


def pair_matching_file_types(mixedlist, ext_a, ext_b):
    """Match pairs of file types with differing extensions."""
    file_dict = {}
    type_a_stems = [str(Path(n).stem) for n in mixedlist if n.find(ext_a) > -1]

    type_a_stems.sort(key=len)
    type_b_stems = [str(Path(n).stem) for n in mixedlist if n.find(ext_b) > -1]
    type_b_stems.sort(key=len)
    if len(type_a_stems) != len(type_b_stems):
        logger.error(
            f"Differing number of {ext_a} ({len(type_b_stems)})"
            + f" and {ext_b} files ({len(type_a_stems)})."
        )
        sys.exit(1)
    for type_b in type_b_stems:
        prefix_len = max(
            [len(prefix([type_b, type_a])) for type_a in type_a_stems]
        )
        match_type_a_idx = [
            i
            for i, type_a in enumerate(type_a_stems)
            if len(prefix([type_b, type_a])) == prefix_len
        ][0]
        match_type_a = type_a_stems.pop(match_type_a_idx)
        type_b_path = [
            Path(p) for p in mixedlist if p.endswith(type_b + "." + ext_b)
        ][0]
        type_a_path = [
            Path(p)
            for p in mixedlist
            if p.endswith(match_type_a + "." + ext_a)
        ][0]
        stem = prefix([type_b, match_type_a])
        file_dict[stem] = {ext_a: type_a_path, ext_b: type_b_path}
    return file_dict


@cli.command()
@click_loguru.init_logger()
@click.option(
    "-s",
    "--shorten_source",
    default=False,
    is_flag=True,
    show_default=True,
    help="Remove invariant dotpaths in source IDs.",
)
@click.argument("setname")
@click.argument("gff_faa_path_list", nargs=-1)
def join_protein_position_info(shorten_source, setname, gff_faa_path_list):
    """Marshal protein and genome sequence information.

    Corresponding GFF and FASTA files must have a corresponding prefix to their
    file names, but theu may occur in any order in the list.  Paths to files
    need not be the same.  Files must be uncompressed. FASTA files must be
    protein files with extension ".faa".  GFF files must have extension ".gff3".

    IDs must correspond between GFF and FASTA files and must be unique across
    the entire set.
    """
    if len(gff_faa_path_list) == 0:
        logger.error("No files in list, exiting.")
        sys.exit(0)
    file_dict = pair_matching_file_types(gff_faa_path_list, GFF_EXT, FAA_EXT)
    set_path = Path(setname)
    set_path.mkdir(parents=True, exist_ok=True)
    file_frame, propfile_path_dict = prepare_protein_files.callback(
        setname, None, stemdict=file_dict
    )
    for stem in file_dict:
        logger.debug(f"Reading GFF file {file_dict[stem][GFF_EXT]}.")
        annotation = gffpd.read_gff3(file_dict[stem][GFF_EXT])
        mrnas = annotation.filter_feature_of_type(
            ["mRNA"]
        ).attributes_to_columns()
        del annotation
        mrnas.drop(
            mrnas.columns.drop(["seq_id", "start", "strand", "ID"]),
            axis=1,
            inplace=True,
        )  # drop non-essential columns
        if shorten_source:
            # drop identical sub-fields in seq_id to keep them visually short (for development)
            split_sources = mrnas["seq_id"].str.split(".", expand=True)
            split_sources = split_sources.drop(
                [
                    i
                    for i in split_sources.columns
                    if len(set(split_sources[i])) == 1
                ],
                axis=1,
            )
            sources = split_sources.agg(".".join, axis=1)
            mrnas["seq_id"] = sources
        mrnas = mrnas.set_index("ID")
        # Make a categorical column, frag_id, based on seq_id
        mrnas["frag_id"] = pd.Categorical(mrnas["seq_id"])
        mrnas.drop(["seq_id"], axis=1, inplace=True)
        propfile_path = propfile_path_dict[stem]
        logger.debug(f"Reading properties file {propfile_path}.")
        prop_frame = pd.read_csv(propfile_path, sep="\t", index_col=0)
        # Drop any mrnas not found in sequence file, e.g., zero-length
        mrnas_len_before = len(mrnas)
        mrnas = mrnas[mrnas.index.isin(prop_frame.index)]
        if len(mrnas) != mrnas_len_before:
            logger.debug(
                f"Dropped {mrnas_len_before - len(mrnas)} "
                + f"mRNA defs from {stem}"
            )
        # sort by largest value
        frag_counts = mrnas["frag_id"].value_counts()
        frag_frame = pd.DataFrame()
        frag_frame["counts"] = frag_counts
        frag_frame["idx"] = range(len(frag_frame))
        frag_frame["frag_id"] = frag_frame.index
        frag_count_path = set_path / f"{stem}-fragment_counts.tsv"
        logger.debug(f"Writing fragment count to file {frag_count_path}")
        frag_frame.set_index(["idx"]).to_csv(frag_count_path, sep="\t")
        del frag_frame
        mrnas["frag_count"] = mrnas["frag_id"].map(frag_counts)
        mrnas.sort_values(
            by=["frag_count", "start"], ascending=[False, True], inplace=True
        )
        frag_id_range = []
        for frag_id in frag_counts.index:
            frag_id_range += list(range(frag_counts[frag_id]))
        mrnas["frag_pos"] = frag_id_range
        del frag_id_range
        mrnas.drop(["frag_count"], axis=1, inplace=True)
        # join GFF info to FASTA info
        joined_path = set_path / f"{stem}-protein_positions.tsv"
        logger.debug(
            f"Joined {len(mrnas)} position and protein info for {stem}"
        )
        mrnas = mrnas.join(prop_frame)
        mrnas.to_csv(joined_path, sep="\t")


@cli.command()
@click_loguru.init_logger()
@click.option(
    "--identity",
    "-i",
    default=0.0,
    help="Minimum sequence ID (0-1). [default: lowest]",
)
@click.option(
    "--clust/--no-clust",
    "-c/-x",
    is_flag=True,
    default=True,
    help="Do cluster calc.",
    show_default=True,
)
@click.argument("setname")
def annotate_homology(identity, clust, setname):
    """Add homology cluster info."""
    set_keys = list(file_frame["stem"])
    concatenated_fasta_name = f"{setname}.faa"
    if clust:
        logger.debug("Doing cluster calculation.")
        cwd = Path.cwd()
        os.chdir(set_path)
        stats, graph, hist, any_, all_ = usearch_cluster.callback(
            concatenated_fasta_name, identity, write_ids=True, delete=False
        )
        os.chdir(cwd)
        del stats, graph, hist, any_, all_
    del fasta_records
    cluster_frame = pd.read_csv(
        set_path / (cluster_set_name(setname, identity) + "-ids.tsv"), sep="\t"
    )
    cluster_frame = cluster_frame.set_index("id")
    logger.debug("Mapping FASTA IDs to cluster properties.")

    def id_to_cluster_property(ident, column):
        try:
            return int(cluster_frame.loc[ident, column])
        except KeyError:
            raise KeyError(f"ID {id} not found in clusters")

    for stem in set_keys:
        frame = frame_dict[stem]
        frame["cluster_id"] = frame.index.map(
            lambda i: id_to_cluster_property(i, "cluster")
        )
        frame["cluster_size"] = frame.index.map(
            lambda i: id_to_cluster_property(i, "siz")
        )
        homology_filename = f"{stem}{HOMOLOGY_ENDING}"
        logger.debug(f"Writing homology file {homology_filename}")
        frame.to_csv(set_path / homology_filename, sep="\t")


@cli.command()
@click_loguru.init_logger()
@click.option("-k", default=6, help="Synteny block length.", show_default=True)
@click.option(
    "-r",
    "--rmer",
    default=False,
    is_flag=True,
    show_default=True,
    help="Allow repeats in block.",
)
@click.argument("setname")
@click.argument("gff_fna_path_list", nargs=-1)
def synteny_anchors(k, rmer, setname, gff_fna_path_list):
    """Calculate synteny anchors."""
    if len(gff_fna_path_list) == 0:
        logger.error("No files in list, exiting.")
        sys.exit(0)
    set_path = Path(setname)
    files_frame, frame_dict = read_files(setname)
    set_keys = list(files_frame["stem"])
    logger.debug(f"Calculating k-mer of length {k} synteny blocks.")
    merge_frame_columns = ["hash", "source"]
    merge_frame = pd.DataFrame(columns=merge_frame_columns)
    for stem in set_keys:
        frame = frame_dict[stem]
        synteny_func_name = synteny_block_func(k, rmer, None, name_only=True)
        frame_len = frame.shape[0]
        map_results = []
        for unused_seq_id, subframe in frame.groupby(by=["seq_id"]):
            hash_closure = synteny_block_func(k, rmer, subframe)
            for i in range(len(subframe)):
                map_results.append(hash_closure(i))
        frame["footprint"] = [map_results[i][0] for i in range(len(frame))]
        frame["hashdir"] = [map_results[i][1] for i in range(len(frame))]
        frame[synteny_func_name] = [
            map_results[i][2] for i in range(len(frame))
        ]
        del map_results
        hash_series = frame[synteny_func_name]
        assigned_hashes = hash_series[hash_series != 0]
        del hash_series
        n_assigned = len(assigned_hashes)
        logger.info(
            f"{stem} has {frame_len} proteins, {n_assigned}"
            + f"of which have {synteny_func_name} hashes,"
        )
        hash_counts = assigned_hashes.value_counts()
        assigned_hash_frame = pd.DataFrame(columns=merge_frame_columns)
        assigned_hash_frame["hash"] = assigned_hashes.unique()
        assigned_hash_frame["source"] = stem
        n_non_unique = n_assigned - len(assigned_hash_frame)
        percent_non_unique = n_non_unique / n_assigned * 100.0
        logger.info(
            f"  of which {n_non_unique} ({percent_non_unique:0.1f})% are non-unique."
        )
        merge_frame.append(assigned_hash_frame)
        del assigned_hash_frame
        # create self_count column in frame
        frame["self_count"] = 0
        for idx, row in frame[frame[synteny_func_name] != 0].iterrows():
            frame.loc[idx, "self_count"] = hash_counts.loc[
                row[synteny_func_name]
            ]
        del hash_counts
    logger.debug(f"Calculating overlap of {len(merge_frame)} hash terms.")
    hash_counts = merge_frame["hash"].value_counts()
    merged_hash_frame = pd.DataFrame(
        index=merge_frame["hash"].unique(), columns=["count"]
    )
    for idx, row in merged_hash_frame.iterrows():
        merged_hash_frame.loc[idx, "count"] = hash_counts.loc[
            row[synteny_func_name]
        ]
    print(f"Merged_hash_frame={merged_hash_frame}")
    merged_hash_frame = merged_hash_frame[merged_hash_frame["count"] > 1]
    print(
        f"after dropping non-matching hashes, len = {len(merged_hash_frame)}"
    )
    print(f"merged hash counts={hash_counts}")
    for stem in set_keys:
        synteny_name = f"{stem}-{synteny_func_name}{SYNTENY_ENDING}"
        logger.debug(
            f"Writing {synteny_func_name} synteny frame {synteny_name}."
        )
        frame_dict[stem].to_csv(set_path / synteny_name, sep="\t")


def dagchainer_id_to_int(ident):
    """Accept DAGchainer ids such as "cl1" and returns an integer."""
    if not ident.startswith("cl"):
        raise ValueError(f"Invalid ID {ident}.")
    id_val = ident[2:]
    if not id_val.isnumeric():
        raise ValueError(f"Non-numeric ID value in {ident}.")
    return int(id_val)


@cli.command()
@click_loguru.init_logger()
@click.argument("setname")
def dagchainer_synteny(setname):
    """Read DAGchainer synteny into homology frames.

    IDs must correspond between DAGchainer files and homology blocks.
    Currently does not calculate DAGchainer synteny.
    """

    cluster_path = Path.cwd() / "out_azulejo" / "clusters.tsv"
    if not cluster_path.exists():
        try:
            azulejo_tool = sh.Command("azulejo_tool")
        except sh.CommandNotFound:
            logger.error("azulejo_tool must be installed first.")
            sys.exit(1)
        logger.debug("Running azulejo_tool clean")
        try:
            output = azulejo_tool(["clean"])
        except sh.ErrorReturnCode:
            logger.error("Error in clean.")
            sys.exit(1)
        logger.debug("Running azulejo_tool run")
        try:
            output = azulejo_tool(["run"])
            print(output)
        except sh.ErrorReturnCode:
            logger.error(
                "Something went wrong in azulejo_tool, check installation."
            )
            sys.exit(1)
        if not cluster_path.exists():
            logger.error(
                "Something went wrong with DAGchainer run.  Please run it manually."
            )
            sys.exit(1)
    synteny_func_name = "dagchainer"
    set_path = Path(setname)
    logger.debug(f"Reading {synteny_func_name} synteny file.")
    synteny_frame = pd.read_csv(
        cluster_path, sep="\t", header=None, names=["cluster", "id"]
    )
    synteny_frame["synteny_id"] = synteny_frame["cluster"].map(
        dagchainer_id_to_int
    )
    synteny_frame = synteny_frame.drop(["cluster"], axis=1)
    cluster_counts = synteny_frame["synteny_id"].value_counts()
    synteny_frame["synteny_count"] = synteny_frame["synteny_id"].map(
        cluster_counts
    )
    synteny_frame = synteny_frame.sort_values(
        by=["synteny_count", "synteny_id"]
    )
    synteny_frame = synteny_frame.set_index(["id"])
    files_frame, frame_dict = read_files(setname)
    set_keys = list(files_frame["stem"])

    def id_to_synteny_property(ident, column):
        try:
            return int(synteny_frame.loc[ident, column])
        except KeyError:
            return 0

    for stem in set_keys:
        homology_frame = frame_dict[stem]
        homology_frame["synteny_id"] = homology_frame.index.map(
            lambda x: id_to_synteny_property(x, "synteny_id")
        )
        homology_frame["synteny_count"] = homology_frame.index.map(
            lambda x: id_to_synteny_property(x, "synteny_count")
        )
        synteny_name = f"{stem}-{synteny_func_name}{SYNTENY_ENDING}"
        logger.debug(
            f"Writing {synteny_func_name} synteny frame {synteny_name}."
        )
        homology_frame.to_csv(set_path / synteny_name, sep="\t")


class ProxySelector:

    """Provide methods for downselection of proxy genes."""

    def __init__(self, frame, prefs):
        """Calculate any joint statistics from frame."""
        self.frame = frame
        self.prefs = prefs
        self.reasons = []
        self.drop_ids = []
        self.first_choice = prefs[0]
        self.first_choice_hits = 0
        self.first_choice_unavailable = 0
        self.cluster_count = 0

    def choose(self, chosen_one, cluster, reason, drop_non_chosen=True):
        """Make the choice, recording stats."""
        self.frame.loc[chosen_one, "reason"] = reason
        self.first_choice_unavailable += int(
            self.first_choice not in set(cluster["stem"])
        )
        self.first_choice_hits += int(
            cluster.loc[chosen_one, "stem"] == self.first_choice
        )
        non_chosen_ones = list(cluster.index)
        non_chosen_ones.remove(chosen_one)
        if drop_non_chosen:
            self.drop_ids += non_chosen_ones
        else:
            self.cluster_count += len(non_chosen_ones)

    def choose_by_preference(
        self, subcluster, cluster, reason, drop_non_chosen=True
    ):
        """Choose in order of preference."""
        stems = subcluster["stem"]
        pref_idxs = [subcluster[stems == pref].index for pref in self.prefs]
        pref_lens = np.array([int(len(idx) > 0) for idx in pref_idxs])
        best_choice = np.argmax(pref_lens)  # first occurrance
        if pref_lens[best_choice] > 1:
            raise ValueError(
                f"subcluster {subcluster} is not unique w.r.t. genome {list(stems)[best_choice]}."
            )
        self.choose(
            pref_idxs[best_choice][0], cluster, reason, drop_non_chosen
        )

    def choose_by_length(self, subcluster, cluster, drop_non_chosen=True):
        """Return an index corresponding to the selected modal/median length."""
        counts = subcluster["protein_len"].value_counts()
        max_count = max(counts)
        if max_count > 1:  # repeated values exist
            max_vals = list(counts[counts == max(counts)].index)
            modal_cluster = subcluster[
                subcluster["protein_len"].isin(max_vals)
            ]
            self.choose_by_preference(
                modal_cluster,
                cluster,
                f"mode{len(modal_cluster)}",
                drop_non_chosen=drop_non_chosen,
            )
        else:
            lengths = list(subcluster["protein_len"])
            median_vals = [
                statistics.median_low(lengths),
                statistics.median_high(lengths),
            ]
            median_pair = subcluster[
                subcluster["protein_len"].isin(median_vals)
            ]
            self.choose_by_preference(
                median_pair, cluster, "median", drop_non_chosen=drop_non_chosen
            )

    def cluster_selector(self, cluster):
        """Calculate which gene in a homology cluster should be left and why."""
        self.cluster_count += 1
        if len(cluster) == 1:
            self.choose(cluster.index[0], cluster, "singleton")
        else:
            for synteny_id, subcluster in cluster.groupby(by=["synteny_id"]):
                if len(subcluster) > 1:
                    self.choose_by_length(
                        subcluster, cluster, drop_non_chosen=(not synteny_id)
                    )
                else:
                    if subcluster["synteny_id"][0] != 0:
                        self.choose(
                            subcluster.index[0],
                            cluster,
                            "bad_synteny",
                            drop_non_chosen=(not synteny_id),
                        )
                    else:
                        self.choose(
                            subcluster.index[0],
                            cluster,
                            "single",
                            drop_non_chosen=(not synteny_id),
                        )

    def downselect_frame(self):
        """Return a frame with reasons for keeping and non-chosen-ones dropped."""
        drop_pct = len(self.drop_ids) * 100.0 / len(self.frame)
        logger.info(
            f"Dropping {len(self.drop_ids)} ({drop_pct:0.1f}%) of {len(self.frame)} genes."
        )
        return self.frame.drop(self.drop_ids)

    def selection_stats(self):
        """Return selection stats."""
        return (
            self.cluster_count,
            self.first_choice_unavailable,
            self.first_choice_hits,
        )


@cli.command()
@click_loguru.init_logger()
@click.argument("setname")
@click.argument("synteny_type")
@click.argument("prefs", nargs=-1)
def proxy_genes(setname, synteny_type, prefs):
    """Calculate a set of proxy genes from synteny files.

    prefs is an optional list of genome stems in order of preference in the proxy calc.
    """
    set_path = Path(setname)
    files_frame, frame_dict = read_files(setname, synteny=synteny_type)
    set_keys = list(files_frame["stem"])
    default_prefs = set_keys.copy()
    default_prefs.reverse()
    if prefs != ():
        for stem in prefs:
            if stem not in default_prefs:
                logger.error(f"Preference {stem} not in {default_prefs}")
                sys.exit(1)
            else:
                default_prefs.remove(stem)
        prefs = list(prefs) + default_prefs
        order = "non-default"
    else:
        prefs = default_prefs
        order = "default"
    logger.debug(
        f"Genome preference for proxy selection in {order} order: {prefs}"
    )
    proxy_frame = None
    for stem in set_keys:
        logger.debug(f"Reading {stem}")
        frame_dict[stem]["stem"] = stem
        if proxy_frame is None:
            proxy_frame = frame_dict[stem]
        else:
            proxy_frame = proxy_frame.append(frame_dict[stem])
    del files_frame
    proxy_frame = proxy_frame.sort_values(
        by=["cluster_size", "cluster_id", "synteny_count", "synteny_id"]
    )
    proxy_filename = f"{setname}-{synteny_type}{PROXY_ENDING}"
    logger.debug(f"Writing initial proxy file {proxy_filename}.")
    proxy_frame.to_csv(set_path / proxy_filename, sep="\t")
    proxy_frame["reason"] = ""
    logger.debug("Downselecting homology clusters.")
    downselector = ProxySelector(proxy_frame, prefs)
    for unused_cluster_id, homology_cluster in proxy_frame.groupby(
        by=["cluster_id"]
    ):  # pylint: disable=unused-variable
        downselector.cluster_selector(homology_cluster)
    downselected = downselector.downselect_frame()
    downselected_filename = (
        f"{setname}-{synteny_type}-downselected{PROXY_ENDING}"
    )
    logger.debug(f"Writing downselected proxy file {downselected_filename}.")
    downselected.to_csv(set_path / downselected_filename, sep="\t")
    # print out stats
    (
        cluster_count,
        first_choice_unavailable,
        first_choice_hits,
    ) = downselector.selection_stats()
    first_choice_percent = (
        first_choice_hits * 100.0 / (cluster_count - first_choice_unavailable)
    )
    first_choice_unavailable_percent = (
        first_choice_unavailable * 100.0 / cluster_count
    )
    logger.info(
        f"First-choice ({prefs[0]}) selections from {cluster_count} homology clusters:"
    )
    logger.info(
        f"   not in cluster: {first_choice_unavailable} ({first_choice_unavailable_percent:.1f}%)"
    )
    logger.info(
        f"   chosen as proxy: {first_choice_hits} ({first_choice_percent:.1f}%)"
    )
