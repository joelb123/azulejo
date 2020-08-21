# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import array
import sys

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import attr
import click
import dask.bag as db
import networkx as nx
import numpy as np
import pandas as pd
import sh
import xxhash
from dask.diagnostics import ProgressBar
from itertools import combinations
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import ANCHOR_HIST_FILE
from .common import CLUSTERS_FILE
from .common import DIRECTIONAL_CATEGORY
from .common import HOMOLOGY_FILE
from .common import PROTEOMOLOGY_FILE
from .common import PROTEOSYN_FILE
from .common import SYNTENY_FILE
from .common import dotpath_to_path
from .common import log_and_add_to_stats
from .common import read_tsv_or_parquet
from .common import remove_tmp_columns
from .common import write_tsv_or_parquet
from .mailboxes import DataMailboxes
from .mailboxes import ExternalMerge

# global constants
CLUSTER_COLS = [
    "syn.anchor_id",
]
MERGE_COLS = ["syn.hash.self_count", "frag.idx"]
DEFAULT_K = 3


@attr.s
class SyntenyBlockHasher(object):

    """Synteny-block hashes via reversible-peatmer method."""

    k = attr.ib(default=5)
    peatmer = attr.ib(default=True)
    prefix = attr.ib(default="syn")

    def hash_name(self, no_prefix=False):
        """Return the string name of the hash function."""
        if no_prefix:
            prefix_str = ""
        else:
            prefix_str = self.prefix + "."
        if self.peatmer:
            return f"{prefix_str}hash.peatmer{self.k}"
        return f"{prefix_str}hash.kmer{self.k}"

    def _hash_kmer(self, kmer):
        """Return a hash of a kmer array."""
        return xxhash.xxh32_intdigest(kmer.tobytes())

    def shingle(self, cluster_series, base, direction, hash):
        """Return a vector of anchor ID's. """
        vec = cluster_series.to_numpy().astype(int)
        steps = np.insert((vec[1:] != vec[:-1]).astype(int), 0, 0).cumsum()
        try:
            assert max(steps) == self.k - 1
        except AssertionError:
            logger.error(
                f"Working around minor error in shingling hash {hash}, base"
                f" {base};"
            )
            logger.error(f"input homology string={vec}")
            logger.error(f"max index = {max(steps)}, should be {self.k-1}")
            steps[np.where(steps > self.k - 1)] = self.k - 1
        if direction == "+":
            return base + steps
        return base + self.k - 1 - steps

    def calculate(self, cluster_series):
        """Return an array of synteny block hashes data."""
        # Maybe the best code I've ever written--JB
        vec = cluster_series.to_numpy().astype(int)
        if self.peatmer:
            uneq_idxs = np.append(np.where(vec[1:] != vec[:-1]), len(vec) - 1)
            runlengths = np.diff(np.append(-1, uneq_idxs))
            positions = np.cumsum(np.append(0, runlengths))[:-1]
            n_mers = len(positions) - self.k + 1
            footprints = pd.array(
                [runlengths[i : i + self.k].sum() for i in range(n_mers)],
                dtype=pd.UInt32Dtype(),
            )
        else:
            n_elements = len(cluster_series)
            n_mers = n_elements - self.k + 1
            positions = np.arange(n_elements)
            footprints = pd.array([self.k] * (n_mers), dtype=pd.UInt32Dtype())
        if n_mers < 1:
            return None
        # Calculate k-mers over indirect index
        kmer_mat = np.array(
            [vec[positions[i : i + self.k]] for i in range(n_mers)]
        )
        fwd_rev_hashes = np.array(
            [
                np.apply_along_axis(self._hash_kmer, 1, kmer_mat),
                np.apply_along_axis(
                    self._hash_kmer, 1, np.flip(kmer_mat, axis=1)
                ),
            ]
        )
        plus_minus = np.array([["+"] * n_mers, ["-"] * n_mers])
        directions = np.take_along_axis(
            plus_minus,
            np.expand_dims(fwd_rev_hashes.argmin(axis=0), axis=0),
            axis=0,
        )[0]
        return pd.DataFrame(
            [
                pd.Categorical(directions, dtype=DIRECTIONAL_CATEGORY),
                footprints,
                pd.array(
                    np.amin(fwd_rev_hashes, axis=0), dtype=pd.UInt32Dtype()
                ),
            ],
            columns=[
                "syn.hash.direction",
                "syn.hash.footprint",
                self.hash_name(),
            ],
            index=cluster_series.index[positions[:n_mers]],
        )


class AmbiguousMerger(object):

    """Counts instance of merges."""

    def __init__(
        self,
        graph_path=Path("synteny.gml"),
        count_key="value",
        ordinal_key="count",
        ambig_key="ambig",
        start_base=0,
        k=None,
        graph=None,
    ):
        """Create arrays as instance attributes."""
        self.graph_path = graph_path
        self.count_key = count_key
        self.ordinal_key = ordinal_key
        self.start_base = start_base
        self.ambig_key = ambig_key
        self.values = array.array("L")
        self.counts = array.array("h")
        self.ambig = array.array("h")
        if graph is None:
            self.graph = nx.Graph()
        else:
            self.graph = graph
        self.k = k

    def _unpack_payloads(self, vec):
        """Unpack TSV ints in payload"""
        wheres = np.where(~vec.mask)[0]
        values = np.array(
            [[int(i) for i in s.split("\t")] for s in vec.compressed()]
        ).transpose()
        return wheres, values

    def merge_func(self, value, count, payload_vec):
        """Return list of merged values."""
        self.values.append(value)
        self.counts.append(count)
        wheres, arr = self._unpack_payloads(payload_vec)
        max_ambig = arr[0].max()
        self.ambig.append(max_ambig)
        if max_ambig == 1:
            self._adjacency_to_graph([f"{i}" for i in arr[1]], value)

    def _adjacency_to_graph(self, nodes, edgename):
        """Turn adjacency data into GML graph file."""
        self.graph.add_nodes_from(nodes)
        edges = combinations(nodes, 2)
        self.graph.add_edges_from(edges, label=f"{edgename}")

    def results(self):
        merge_frame = pd.DataFrame(
            {self.count_key: self.counts, self.ambig_key: self.ambig},
            index=self.values,
            dtype=pd.UInt32Dtype(),
        )
        merge_frame.sort_values(
            by=[self.ambig_key, self.count_key], inplace=True
        )
        unambig_frame = merge_frame[merge_frame[self.ambig_key] == 1].copy()
        n_unambig = len(unambig_frame)
        unambig_frame[self.ordinal_key] = (
            pd.array(
                range(self.start_base, self.start_base + n_unambig),
                dtype=pd.UInt32Dtype(),
            )
            * self.k
        )
        del unambig_frame[self.ambig_key]

        ambig_frame = merge_frame[merge_frame[self.ambig_key] > 1].copy()
        ambig_frame = ambig_frame.rename(
            columns={self.count_key: self.count_key + ".ambig"}
        )
        ambig_frame[self.ordinal_key + ".ambig"] = (
            pd.array(
                range(
                    self.start_base + n_unambig,
                    self.start_base + len(ambig_frame) + n_unambig,
                ),
                dtype=pd.UInt32Dtype(),
            )
            * self.k
        )
        del merge_frame
        nx.write_gml(self.graph, self.graph_path)
        return unambig_frame, ambig_frame, self.graph


@cli.command()
@click_loguru.init_logger()
@click.option(
    "-k", default=DEFAULT_K, help="Synteny block length.", show_default=True
)
@click.option(
    "--peatmer/--kmer",
    default=True,
    is_flag=True,
    show_default=True,
    help="Allow repeats in block.",
)
@click.option(
    "--greedy/--no-greedy",
    is_flag=True,
    default=True,
    show_default=True,
    help="Assign ambiguous hashes to longest frag.",
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    default=True,
    show_default=True,
    help="Process in parallel.",
)
@click.argument("setname")
def synteny_anchors(k, peatmer, setname, parallel, greedy):
    """Calculate synteny anchors."""
    if k < 2:
        logger.error("k must be at least 2.")
        sys.exit(1)
    options = click_loguru.get_global_options()
    set_path = Path(setname)
    file_stats_path = set_path / PROTEOMOLOGY_FILE
    proteomes = read_tsv_or_parquet(file_stats_path)
    n_proteomes = len(proteomes)
    clusters = read_tsv_or_parquet(set_path / CLUSTERS_FILE)
    n_clusters = len(clusters)
    hasher = SyntenyBlockHasher(k=k, peatmer=peatmer)
    hash_mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "hash_merge"),
    )
    hash_mb.write_headers("hash\n")
    arg_list = []
    hash_stats_list = []
    for idx, row in proteomes.iterrows():
        arg_list.append((idx, row["path"],))
    if parallel:
        bag = db.from_sequence(arg_list)
    if not options.quiet:
        logger.info(
            f"Calculating {hasher.hash_name(no_prefix=True)} synteny anchors"
            + f" for {n_proteomes} proteomes"
        )
        ProgressBar().register()
    if parallel:
        hash_stats_list = bag.map(
            calculate_synteny_hashes, mailboxes=hash_mb, hasher=hasher
        ).compute()
    else:
        for args in arg_list:
            hash_stats_list.append(
                calculate_synteny_hashes(
                    args, mailboxes=hash_mb, hasher=hasher
                )
            )
    hash_stats = (
        pd.DataFrame.from_dict(hash_stats_list).set_index("idx").sort_index()
    )
    proteomes = log_and_add_to_stats(proteomes, hash_stats)
    del hash_stats_list, hash_stats
    merger = ExternalMerge(
        file_path_func=hash_mb.path_to_mailbox, n_merge=n_proteomes
    )
    merger.init("hash")
    merge_counter = AmbiguousMerger(
        count_key="syn.hash.ortho_count",
        ordinal_key="tmp.base",
        ambig_key="tmp.max_ambig",
        graph_path=set_path / "synteny.gml",
        k=k,
    )
    unambig_hashes, ambig_hashes, fragment_synteny_graph = merger.merge(
        merge_counter
    )
    last_anchor_base = len(unambig_hashes) + len(ambig_hashes) + 1
    hash_mb.delete()
    del merger, merge_counter
    ambig_mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "ambig_merge"),
    )
    ambig_mb.write_headers("hash\n")
    synteny_stats_list = []
    if not options.quiet:
        logger.info(
            f"Merging {len(unambig_hashes)} unambiguous "
            + f"and disambiguating {len(ambig_hashes)} synteny anchors into "
            + f"{n_proteomes} proteomes"
        )
        ProgressBar().register()
    if parallel:
        synteny_stats_list = bag.map(
            merge_unambig_hashes,
            unambig_hashes=unambig_hashes,
            ambig_hashes=ambig_hashes,
            hasher=hasher,
            ambig_mb=ambig_mb,
        ).compute()
    else:
        for args in arg_list:
            synteny_stats_list.append(
                merge_unambig_hashes(
                    args,
                    unambig_hashes=unambig_hashes,
                    ambig_hashes=ambig_hashes,
                    hasher=hasher,
                    ambig_mb=ambig_mb,
                )
            )
    synteny_stats = (
        pd.DataFrame.from_dict(synteny_stats_list)
        .set_index("idx")
        .sort_index()
    )
    proteomes = log_and_add_to_stats(proteomes, synteny_stats)
    del synteny_stats_list, synteny_stats, unambig_hashes
    #
    logger.info(
        f"Reducing {proteomes['tmp.hashes.pending_disambig'].max()} disambiguation hashes"
        + " via external merge"
    )
    disambig_merger = ExternalMerge(
        file_path_func=ambig_mb.path_to_mailbox, n_merge=n_proteomes
    )
    disambig_merger.init("hash")
    disambig_merge_counter = AmbiguousMerger(
        count_key="syn.hash.disambig.ortho_count",
        ordinal_key="tmp.disambig.base",
        ambig_key="tmp.disambig.max_ambig",
        graph_path=set_path / "synteny-disambig.gml",
        start_base=last_anchor_base,
        k=k,
        graph=fragment_synteny_graph,
    )
    (
        disambig_hashes,
        still_ambig_hashes,
        fragment_synteny_graph,
    ) = disambig_merger.merge(disambig_merge_counter)
    ambig_mb.delete()
    #
    cluster_mb = DataMailboxes(
        n_boxes=n_clusters,
        mb_dir_path=(set_path / "mailboxes" / "anchor_merge"),
        file_extension="tsv",
    )
    cluster_mb.write_tsv_headers(CLUSTER_COLS)
    if parallel:
        bag = db.from_sequence(arg_list)
    disambig_stats_list = []
    if not options.quiet:
        if greedy:
            greedy_txt = "Greedy-merging"
        else:
            greedy_txt = "Merging"
        logger.info(
            f"{greedy_txt} {len(disambig_hashes)} disambiguated "
            + f"and {len(still_ambig_hashes)} partly-ambiguous synteny anchors into {n_proteomes}"
            " proteomes"
        )
        ProgressBar().register()
    if parallel:
        disambig_stats_list = bag.map(
            merge_disambig_hashes,
            disambig_hashes=disambig_hashes,
            still_ambig_hashes=still_ambig_hashes,
            ambig_hashes=ambig_hashes,
            hasher=hasher,
            n_proteomes=n_proteomes,
            cluster_writer=cluster_mb.locked_open_for_write,
            greedy=greedy,
        ).compute()
    else:
        for args in arg_list:
            disambig_stats_list.append(
                merge_disambig_hashes(
                    args,
                    disambig_hashes=disambig_hashes,
                    still_ambig_hashes=still_ambig_hashes,
                    ambig_hashes=ambig_hashes,
                    hasher=hasher,
                    n_proteomes=n_proteomes,
                    cluster_writer=cluster_mb.locked_open_for_write,
                    greedy=greedy,
                )
            )
    disambig_stats = (
        pd.DataFrame.from_dict(disambig_stats_list)
        .set_index("idx")
        .sort_index()
    )
    proteomes = log_and_add_to_stats(proteomes, disambig_stats)
    del disambig_stats_list, disambig_stats
    write_tsv_or_parquet(proteomes, set_path / PROTEOSYN_FILE)
    # merge anchor info into clusters
    arg_list = [(i,) for i in range(n_clusters)]
    if parallel:
        bag = db.from_sequence(arg_list)
    else:
        anchor_stats = []
    if not options.quiet:
        logger.info(f"Joining anchor info to {n_clusters} clusters:")
        ProgressBar().register()
    if parallel:
        anchor_stats = bag.map(
            join_synteny_to_clusters,
            mailbox_reader=cluster_mb.open_then_delete,
            cluster_parent=set_path / "homology",
        ).compute()
    else:
        for args in arg_list:
            anchor_stats.append(
                join_synteny_to_clusters(
                    args,
                    mailbox_reader=cluster_mb.open_then_delete,
                    cluster_parent=set_path / "homology",
                )
            )
    cluster_mb.delete()
    anchor_frame = pd.DataFrame.from_dict(anchor_stats)
    anchor_frame.set_index(["clust_id"], inplace=True)
    anchor_frame.sort_index(inplace=True)
    # with pd.option_context(
    #    "display.max_rows", None, "display.float_format", "{:,.2f}%".format
    # ):
    #    logger.info(anchor_frame)
    proteomes = concat_without_overlap(clusters, anchor_frame)
    write_tsv_or_parquet(
        proteomes, set_path / CLUSTERS_FILE, float_format="%5.2f"
    )
    mean_gene_synteny = (
        proteomes["in_synteny"].sum() * 100.0 / proteomes["size"].sum()
    )
    mean_clust_synteny = proteomes["synteny_pct"].mean()
    logger.info(
        f"Mean anchor coverage: {mean_gene_synteny: .1f}% (on proteins)"
    )
    logger.info(
        f"Mean cluster anchor coverage: {mean_clust_synteny:.1f}% (on clusters)"
    )


def concat_without_overlap(df1, df2):
    """Concatenate two frames on columns, deleting any overlapping columns first."""
    overlapping = set(df1.columns).intersection(df2.columns)
    if len(overlapping) > 0:
        df1 = df1.drop(columns=overlapping)
    return pd.concat([df1, df2], axis=1)


def calculate_synteny_hashes(args, mailboxes=None, hasher=None):
    """Calculate synteny hashes for proteins per-fragment."""
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    hom = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
    hom["tmp.nan_group"] = (
        (hom["hom.cluster"].isnull()).astype(int).cumsum() + 1
    ) * (~hom["hom.cluster"].isnull())
    hom.replace(to_replace={"tmp.nan_group": 0}, value=pd.NA, inplace=True)
    hash_name = hasher.hash_name()
    syn_list = []
    for unused_id_tuple, subframe in hom.groupby(
        by=["frag.id", "tmp.nan_group"]
    ):
        syn_list.append(hasher.calculate(subframe["hom.cluster"]))
    syn = hom.join(
        pd.concat([df for df in syn_list if df is not None], axis=0)
    )
    del syn_list
    syn["tmp.i"] = pd.array(range(len(syn)), dtype=pd.UInt32Dtype())
    hash_counts = syn[hash_name].value_counts()
    syn["syn.hash.self_count"] = pd.array(
        syn[hash_name].map(hash_counts), dtype=pd.UInt32Dtype()
    )
    hash_frag_count_arr = pd.array([pd.NA] * len(syn), dtype=pd.UInt32Dtype())
    hash_is_null = syn[hash_name].isnull()
    for unused_frag, subframe in syn.groupby(by=["frag.id"]):
        try:
            frag_hash_counts = subframe[hash_name].value_counts()
        except ValueError:
            continue
        for unused_i, row in subframe.iterrows():
            row_no = row["tmp.i"]
            if not hash_is_null[row_no]:
                hash_val = row[hash_name]
                hash_frag_count_arr[row_no] = frag_hash_counts[hash_val]
    syn["syn.hash.frag_count"] = hash_frag_count_arr
    del syn["tmp.i"]
    write_tsv_or_parquet(syn, outpath / SYNTENY_FILE, remove_tmp=False)
    unique_hashes = (
        syn[[hash_name] + MERGE_COLS]
        .drop_duplicates(subset=[hash_name])
        .dropna(how="any")
    )
    unique_hashes = unique_hashes.set_index(hash_name).sort_index()
    with mailboxes.locked_open_for_write(idx) as fh:
        unique_hashes.to_csv(fh, header=False, sep="\t")
    in_hash = syn[hash_name].notna().sum()
    n_assigned = syn["hom.cluster"].notna().sum()
    hash_pct = in_hash * 100.0 / n_assigned
    hash_stats = {
        "idx": idx,
        "path": dotpath,
        "hom.clusters": n_assigned,
        "syn.hashes": in_hash,
        "syn.hash_pct": hash_pct,
    }
    return hash_stats


def merge_unambig_hashes(
    args, unambig_hashes=None, ambig_hashes=None, hasher=None, ambig_mb=None,
):
    """Merge unambiguous synteny hashes into proteomes per-proteome."""
    hash_name = hasher.hash_name()
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = join_on_col_with_NA(syn, unambig_hashes, hash_name)
    syn = join_on_col_with_NA(syn, ambig_hashes, hash_name)
    n_proteins = len(syn)
    syn["tmp.i"] = pd.array(range(len(syn)), dtype=pd.UInt32Dtype())
    anchor_blocks = np.array([np.nan] * n_proteins)
    for hash_val, subframe in syn.groupby(by=["tmp.base"]):
        # Note that base values are ordered with lower ortho counts first
        for unused_i, row in subframe.iterrows():
            footprint = row["syn.hash.footprint"]
            row_no = row["tmp.i"]
            anchor_blocks[row_no : row_no + footprint] = hasher.shingle(
                syn["hom.cluster"][row_no : row_no + footprint],
                row["tmp.base"],
                row["syn.hash.direction"],
                row[hash_name],
            )
    syn["syn.anchor_id"] = pd.array(anchor_blocks, dtype=pd.UInt32Dtype())
    # Calculate disambiguation hashes and write them out for merge
    disambig_frame_list = []
    for unused_frag, subframe in syn.groupby(by=["frag.id"]):
        disambig_frame_list.append(calculate_disambig_hashes(subframe))
    disambig_fr = pd.concat(
        [df for df in disambig_frame_list if df is not None]
    )
    disambig_fr = disambig_fr.dropna(how="all")
    syn = syn.join(disambig_fr)
    write_tsv_or_parquet(syn, outpath / SYNTENY_FILE, remove_tmp=False)
    upstream_hashes = (
        syn[["syn.disambig_upstr"] + MERGE_COLS]
        .dropna(how="any")
        .rename(columns={"syn.disambig_upstr": "hash"})
    )
    downstream_hashes = (
        syn[["syn.disambig_downstr"] + MERGE_COLS]
        .dropna(how="any")
        .rename(columns={"syn.disambig_downstr": "hash"})
    )
    unique_hashes = pd.concat(
        [upstream_hashes, downstream_hashes], ignore_index=True
    ).drop_duplicates(subset=["hash"])
    unique_hashes = unique_hashes.set_index("hash").sort_index()
    with ambig_mb.locked_open_for_write(idx) as fh:
        unique_hashes.to_csv(fh, header=False, sep="\t")
    # Do some synteny stats
    in_hash = syn[hash_name].notna().sum()
    n_unambig = syn["tmp.base"].notna().sum()
    n_ambig = syn["tmp.base.ambig"].notna().sum()
    unambig_pct = n_unambig * 100.0 / in_hash
    ambig_pct = n_ambig * 100.0 / in_hash
    in_synteny = syn["syn.anchor_id"].notna().sum()
    synteny_pct = in_synteny * 100.0 / in_synteny
    synteny_stats = {
        "idx": idx,
        "path": dotpath,
        "syn.hashes.n": in_hash,
        "syn.hashes.unambig": n_unambig,
        "syn.hashes.unambig_pct": unambig_pct,
        "syn.hashes.ambig": n_ambig,
        "syn.hashes.ambig_pct": ambig_pct,
        "tmp.prot.in_synteny": in_synteny,
        "tmp.hashes.pending_disambig": len(disambig_fr),
    }
    return synteny_stats


def join_on_col_with_NA(left, right, col_name):
    """Join on a temporary column of type 'O'."""
    tmp_col_name = "tmp." + col_name
    left[tmp_col_name] = left[col_name].astype("O")
    merged = pd.merge(
        left, right, left_on=tmp_col_name, right_index=True, how="left"
    )
    del merged[tmp_col_name]
    return merged


def merge_disambig_hashes(
    args,
    disambig_hashes=None,
    still_ambig_hashes=None,
    ambig_hashes=None,
    hasher=None,
    n_proteomes=None,
    cluster_writer=None,
    greedy=False,
):
    """Merge disambiguated synteny hashes into proteomes per-proteome."""
    plain_hash_name = hasher.hash_name(no_prefix=True)
    hash_name = "syn." + plain_hash_name
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = join_on_col_with_NA(syn, disambig_hashes, "syn.disambig_upstr")
    syn = join_on_col_with_NA(syn, disambig_hashes, "syn.disambig_downstr")
    syn = join_on_col_with_NA(syn, still_ambig_hashes, "syn.disambig_upstr")
    syn = join_on_col_with_NA(syn, still_ambig_hashes, "syn.disambig_downstr")
    del syn["syn.disambig_upstr"], syn["syn.disambig_downstr"]
    for dup_col in [
        "syn.hash.disambig.ortho_count",
        "syn.hash.disambig.ortho_count.ambig",
        "tmp.disambig.base",
        "tmp.disambig.base.ambig",
        "tmp.disambig.max_ambig",
    ]:
        xcol = dup_col + "_x"
        ycol = dup_col + "_y"
        syn[dup_col] = syn[xcol].fillna(syn[ycol])
        del syn[xcol], syn[ycol]
    syn["syn.hash.ortho_count"] = syn["syn.hash.ortho_count"].fillna(
        syn["syn.hash.disambig.ortho_count"]
    )
    del syn["syn.hash.disambig.ortho_count"]
    n_proteins = len(syn)
    syn["tmp.i"] = range(len(syn))
    # Do disambiguated fills
    disambig_fills = np.array([np.nan] * n_proteins)
    for hash_val, subframe in syn.groupby(by=["tmp.disambig.base"]):
        for unused_i, row in subframe.iterrows():
            footprint = row["syn.hash.footprint"]
            row_no = row["tmp.i"]
            disambig_fills[row_no : row_no + footprint] = hasher.shingle(
                syn["hom.cluster"][row_no : row_no + footprint],
                row["tmp.disambig.base"],
                row["syn.hash.direction"],
                row[hash_name],
            )
    syn["syn.anchor_id"] = syn["syn.anchor_id"].fillna(
        pd.Series(disambig_fills, index=syn.index)
    )
    n_disambiguated = (syn["syn.anchor_id"] == disambig_fills).sum()
    # Deal with ambiguous hashes by adding non-ambiguous examples
    nonambig_fills = np.array([np.nan] * n_proteins)
    for hash_val, subframe in syn.groupby(by=["tmp.base.ambig"]):
        if not greedy and len(subframe) > 1:
            continue
        for unused_i, row in subframe.iterrows():
            footprint = row["syn.hash.footprint"]
            row_no = row["tmp.i"]
            nonambig_fills[row_no : row_no + footprint] = hasher.shingle(
                syn["hom.cluster"][row_no : row_no + footprint],
                row["tmp.base.ambig"],
                row["syn.hash.direction"],
                row[hash_name],
            )
            break
    syn["syn.anchor_id"] = syn["syn.anchor_id"].fillna(
        pd.Series(nonambig_fills, index=syn.index)
    )
    n_nonambig = (syn["syn.anchor_id"] == nonambig_fills).sum()
    # Delete temporary columns
    non_needed_cols = [
        "tmp.i",
        "syn.hash.disambig.ortho_count.ambig",
        "syn.hash.ortho_count.ambig",
        "syn.hash.frag_count",
        "syn.hash.footprint",
        "syn.hash.direction",
        "syn.hash.self_count",
    ]
    syn = syn.drop(columns=non_needed_cols)
    write_tsv_or_parquet(
        syn, outpath / SYNTENY_FILE,
    )
    for cluster_id, subframe in syn.groupby(by=["hom.cluster"]):
        with cluster_writer(cluster_id) as fh:
            subframe[CLUSTER_COLS].dropna().to_csv(fh, header=False, sep="\t")
    n_assigned = n_proteins - syn["hom.cluster"].isnull().sum()
    # Do histogram of blocks
    anchor_counts = syn["syn.anchor_id"].value_counts()
    anchor_hist = pd.DataFrame(anchor_counts.value_counts()).sort_index()
    anchor_hist = anchor_hist.rename(
        columns={"syn.anchor_id": "hash.self_count"}
    )
    anchor_hist["pct_anchors"] = (
        anchor_hist["hash.self_count"] * anchor_hist.index * 100.0 / n_assigned
    )
    write_tsv_or_parquet(anchor_hist, outpath / ANCHOR_HIST_FILE)
    # Do histogram of anchors
    in_synteny = syn["syn.anchor_id"].notna().sum()
    n_assigned = syn["hom.cluster"].notna().sum()
    avg_ortho = syn["syn.hash.ortho_count"].mean()
    synteny_pct = in_synteny * 100.0 / n_assigned
    unassigned_pct = (n_assigned - in_synteny) * 100.0 / n_assigned
    synteny_stats = {
        "idx": idx,
        "path": dotpath,
        "hom.clusters": n_assigned,
        "syn.anchors": in_synteny,
        "syn.pct": synteny_pct,
        "syn.disambig": n_disambiguated,
        "syn.nonambig": n_nonambig,
        "syn.uassigned": n_assigned - in_synteny,
        "syn.unassigned_pct": unassigned_pct,
        "syn.fom": avg_ortho * 100.0 / n_proteomes,
    }
    return synteny_stats


def calculate_disambig_hashes(df):
    """Calculate disambiguation hashes per-fragment."""
    hash2_fr = df[["syn.anchor_id", "tmp.base.ambig"]].copy()
    hash2_fr = hash2_fr.rename(columns={"syn.anchor_id": "tmp.anchor_id"})
    hash2_fr["tmp.upstream"] = fill_na_with_last_valid(df["syn.anchor_id"])
    hash2_fr["tmp.downstream"] = fill_na_with_last_valid(
        df["syn.anchor_id"], flip=True
    )
    hash2_fr["tmp.i"] = range(len(hash2_fr))
    upstream_hash = pd.array([pd.NA] * len(hash2_fr), dtype=pd.UInt32Dtype())
    downstream_hash = pd.array([pd.NA] * len(hash2_fr), dtype=pd.UInt32Dtype())
    hash2_fr["syn.disambig_upstr"] = pd.NA
    hash2_fr["syn.disambig_downstr"] = pd.NA
    for unused_id, row in hash2_fr.iterrows():
        row_no = row["tmp.i"]
        ambig_base = row["tmp.base.ambig"]
        upstream_unambig = row["tmp.upstream"]
        downstream_unambig = row["tmp.downstream"]
        if pd.notna(ambig_base):
            if pd.notna(upstream_unambig):
                upstream_hash[row_no] = xxhash.xxh32_intdigest(
                    np.array([upstream_unambig, ambig_base]).tobytes()
                )
            if pd.notna(downstream_unambig):
                downstream_hash[row_no] = xxhash.xxh32_intdigest(
                    np.array([ambig_base, downstream_unambig]).tobytes()
                )
    hash2_fr["syn.disambig_upstr"] = upstream_hash
    hash2_fr["syn.disambig_downstr"] = downstream_hash
    hash2_fr = remove_tmp_columns(hash2_fr)
    return hash2_fr


def fill_na_with_last_valid(ser, flip=False):
    """Input a series with NA values, returns a series with those values filled."""
    vec = ser.isnull().to_numpy()
    if flip:
        vec = np.flip(vec)
    uneq_idxs = np.append(np.where(vec[1:] != vec[:-1]), len(vec) - 1)
    if len(uneq_idxs) == 0:
        return ser
    runlengths = np.diff(np.append(-1, uneq_idxs))
    positions = np.cumsum(np.append(0, runlengths))[:-1]
    first_null_idxs = np.where(vec[positions])
    first_null_pos = positions[first_null_idxs]
    null_runs = runlengths[first_null_idxs]
    if flip:
        fill_vals = np.append(pd.NA, np.flip(ser.to_numpy()))[first_null_pos]
    else:
        fill_vals = np.append(pd.NA, ser.to_numpy())[first_null_pos]
    out_arr = pd.array([pd.NA] * len(ser), dtype=pd.UInt32Dtype(),)
    for i, pos in enumerate(first_null_pos):
        for j in range(null_runs[i]):
            out_arr[pos + j] = fill_vals[i]
    if flip:
        out_arr = np.flip(out_arr)
    out_ser = pd.Series(out_arr, index=ser.index)
    return out_ser


def join_synteny_to_clusters(args, cluster_parent=None, mailbox_reader=None):
    """Read homology info from mailbox and join it to proteome file."""
    idx = args[0]
    cluster_path = cluster_parent / f"{idx}.parq"
    cluster = pd.read_parquet(cluster_path)
    n_cluster = len(cluster)
    with mailbox_reader(idx) as fh:
        synteny_frame = pd.read_csv(fh, sep="\t", index_col=0).convert_dtypes()
        in_synteny = len(synteny_frame)
    # delete columns from previous merge
    for col in synteny_frame.columns:
        if col in cluster.columns:
            del cluster[col]
    clust_syn = concat_without_overlap(cluster, synteny_frame)
    write_tsv_or_parquet(clust_syn, cluster_path)
    return {
        "clust_id": idx,
        "in_synteny": in_synteny,
        "synteny_pct": in_synteny * 100.0 / n_cluster,
    }


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
                "Something went wrong with DAGchainer run.  Please run it"
                " manually."
            )
            sys.exit(1)
    synteny_hash_name = "dagchainer"
    set_path = Path(setname)
    logger.debug(f"Reading {synteny_hash_name} synteny file.")
    syn = pd.read_csv(
        cluster_path, sep="\t", header=None, names=["hom.cluster", "id"]
    )
    syn["synteny_id"] = syn["hom.cluster"].map(dagchainer_id_to_int)
    syn = syn.drop(["hom.cluster"], axis=1)
    cluster_counts = syn["synteny_id"].value_counts()
    syn["synteny_count"] = syn["synteny_id"].map(cluster_counts)
    syn = syn.sort_values(by=["synteny_count", "synteny_id"])
    syn = syn.set_index(["id"])
    files_frame, frame_dict = read_files(setname)
    set_keys = list(files_frame["stem"])

    def id_to_synteny_property(ident, column):
        try:
            return int(syn.loc[ident, column])
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
        synteny_name = f"{stem}-{synteny_hash_name}{SYNTENY_ENDING}"
        logger.debug(
            f"Writing {synteny_hash_name} synteny frame {synteny_name}."
        )
        homology_frame.to_csv(set_path / synteny_name, sep="\t")
