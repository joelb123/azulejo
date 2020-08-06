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
from .common import read_tsv_or_parquet
from .common import write_tsv_or_parquet
from .mailboxes import DataMailboxes
from .mailboxes import ExternalMerge

# global constants
CLUSTER_COLS = [
    "syn.anchor_id",
    "syn.self_count",
    "syn.frag_count",
    "syn.direction",
    "syn.ortho_count",
]

MERGE_COLS = ["tmp.self_count", "frag.idx"]


@attr.s
class SyntenyBlockHasher(object):

    """Synteny-block hashes via reversible-peatmer method."""

    k = attr.ib(default=5)
    peatmer = attr.ib(default=True)
    prefix = attr.ib(default="tmp")

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
                    np.amin(fwd_rev_hashes, axis=0), dtype=pd.Int64Dtype()
                ),
            ],
            columns=["tmp.direction", "tmp.footprint", self.hash_name()],
            index=cluster_series.index[positions[:n_mers]],
        )


@attr.s
class SimpleMerger(object):

    """Counts instance of merges."""

    count_key = attr.ib(default="value")
    ordinal_key = attr.ib(default="count")
    values = array.array("L")
    counts = array.array("h")

    def merge_func(self, value, count, unused_payload_vec):
        """Return list of merged values."""
        self.values.append(value)
        self.counts.append(count)

    def results(self):
        merge_frame = pd.DataFrame(
            pd.array(self.counts, dtype=pd.UInt32Dtype()),
            columns=[self.count_key],
            index=self.values,
        )
        merge_frame.sort_values(by=[self.count_key], inplace=True)
        merge_frame[self.ordinal_key] = pd.array(
            range(len(merge_frame)), dtype=pd.UInt32Dtype()
        )
        return merge_frame


@attr.s
class AmbiguousMerger(object):

    """Counts instance of merges."""

    graph_path = attr.ib(default=Path("synteny.gml"))
    count_key = attr.ib(default="value")
    ordinal_key = attr.ib(default="count")
    ambig_key = attr.ib(default="ambig")
    values = array.array("L")
    counts = array.array("h")
    ambig = array.array("h")
    graph = nx.Graph()

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
        self.ambig.append(arr[0].max())
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
        unambig_frame[self.ordinal_key] = pd.array(
            range(len(unambig_frame)), dtype=pd.UInt32Dtype()
        )
        del unambig_frame[self.ambig_key]
        ambig_frame = merge_frame[merge_frame[self.ambig_key] > 1].copy()
        ambig_frame = ambig_frame.rename(
            columns={self.count_key: self.count_key + ".ambig"}
        )
        ambig_frame[self.ordinal_key + ".ambig"] = pd.array(
            range(len(ambig_frame)), dtype=pd.UInt32Dtype()
        )
        del merge_frame
        logger.info(f"Writing fragment synteny graph to {self.graph_path}")
        nx.write_gml(self.graph, self.graph_path)
        return unambig_frame, ambig_frame, self.graph


@cli.command()
@click_loguru.init_logger()
@click.option("-k", default=5, help="Synteny block length.", show_default=True)
@click.option(
    "--peatmer/--kmer",
    default=True,
    is_flag=True,
    show_default=True,
    help="Allow repeats in block.",
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    default=True,
    show_default=True,
    help="Process in parallel.",
)
@click.argument("setname")
def synteny_anchors(k, peatmer, setname, parallel):
    """Calculate synteny anchors."""
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
    cluster_mb = DataMailboxes(
        n_boxes=n_clusters,
        mb_dir_path=(set_path / "mailboxes" / "anchor_merge"),
        file_extension="tsv",
    )
    cluster_mb.write_tsv_headers(CLUSTER_COLS)
    arg_list = []
    n_hashes_list = []
    for idx, row in proteomes.iterrows():
        arg_list.append((idx, row["path"],))
    if parallel:
        bag = db.from_sequence(arg_list)
    if not options.quiet:
        logger.info(
            "Calculating synteny anchors using the"
            f" {hasher.hash_name(no_prefix=True)} function"
            + f" for {n_proteomes} proteomes"
        )
        ProgressBar().register()
    if parallel:
        n_hashes_list = bag.map(
            calculate_synteny_hashes, mailboxes=hash_mb, hasher=hasher
        ).compute()
    else:
        for args in arg_list:
            n_hashes_list.append(
                calculate_synteny_hashes(
                    args, mailboxes=hash_mb, hasher=hasher
                )
            )
    logger.info(f"Reducing {sum(n_hashes_list)} hashes via external merge")
    merger = ExternalMerge(
        file_path_func=hash_mb.path_to_mailbox, n_merge=n_proteomes
    )
    merger.init("hash")
    merge_counter = AmbiguousMerger(
        count_key="tmp.ortho_count",
        ordinal_key="tmp.base",
        ambig_key="tmp.max_ambig",
        graph_path=set_path / "synteny.gml",
    )
    unambig_hashes, ambig_hashes, graph = merger.merge(merge_counter)
    unambig_hashes["tmp.base"] *= k
    ambig_hashes["tmp.base.ambig"] *= k
    hash_mb.delete()
    ret_list = []
    if not options.quiet:
        logger.info(
            f"Merging {len(unambig_hashes)} synteny anchors into {n_proteomes}"
            " proteomes"
        )
        ProgressBar().register()
    if parallel:
        ret_list = bag.map(
            merge_synteny_hashes,
            unambig_hashes=unambig_hashes,
            ambig_hashes=ambig_hashes,
            hasher=hasher,
            n_proteomes=n_proteomes,
            file_writer=cluster_mb.locked_open_for_write,
        ).compute()
    else:
        for args in arg_list:
            ret_list.append(
                merge_synteny_hashes(
                    args,
                    unambig_hashes=unambig_hashes,
                    ambig_hashes=ambig_hashes,
                    hasher=hasher,
                    n_proteomes=n_proteomes,
                    file_writer=cluster_mb.locked_open_for_write,
                )
            )
    synteny_stats = pd.DataFrame.from_dict(ret_list)
    synteny_stats = synteny_stats.set_index("idx").sort_index()
    del synteny_stats["syn.unamb_pct"]
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.float_format",
        "{:,.1f}%".format,
    ):
        logger.info(synteny_stats)
    del synteny_stats["path"], synteny_stats["hom.clusters"]
    proteomes = pd.concat([proteomes, synteny_stats], axis=1)
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
    logger.info(f"Mean anchor coverage on genes: {mean_gene_synteny: .1f}%")
    logger.info(f"Mean anchor coverage on clusters: {mean_clust_synteny:.1f}%")


def concat_without_overlap(df1, df2):
    """Concatenate two frames on columns, deleting any overlapping columns first."""
    overlapping = set(df1.columns).intersection(df2.columns)
    if len(overlapping) > 0:
        df1 = df1.drop(columns=overlapping)
    return pd.concat([df1, df2], axis=1)


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


def calculate_synteny_hashes(args, mailboxes=None, hasher=None):
    """Calculate synteny hashes for protein genes."""
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
    syn = pd.concat([df for df in syn_list if df is not None], axis=0)
    del syn_list
    syn["tmp.frag.id"] = syn.index.map(hom["frag.id"])
    syn["tmp.i"] = pd.array(range(len(syn)), dtype=pd.UInt32Dtype())
    hash_counts = syn[hash_name].value_counts()
    syn["tmp.self_count"] = pd.array(
        syn[hash_name].map(hash_counts), dtype=pd.UInt32Dtype()
    )
    frag_count_arr = pd.array([pd.NA] * len(syn), dtype=pd.UInt32Dtype())
    hash_is_null = syn[hash_name].isnull()
    for unused_frag, subframe in syn.groupby(by=["tmp.frag.id"]):
        try:
            frag_hash_counts = subframe[hash_name].value_counts()
        except ValueError:
            continue
        for unused_i, row in subframe.iterrows():
            row_no = row["tmp.i"]
            if not hash_is_null[row_no]:
                hash_val = row[hash_name]
                frag_count_arr[row_no] = frag_hash_counts[hash_val]
    syn["tmp.frag_count"] = frag_count_arr
    del syn["tmp.i"]
    write_tsv_or_parquet(syn, outpath / SYNTENY_FILE, remove_tmp=False)
    syn["frag.idx"] = hom["frag.idx"]
    unique_hashes = syn[[hash_name] + MERGE_COLS].drop_duplicates(
        subset=[hash_name]
    )
    unique_hashes = unique_hashes.set_index(hash_name).sort_index()
    with mailboxes.locked_open_for_write(idx) as fh:
        unique_hashes.to_csv(fh, header=False, sep="\t")
    return len(unique_hashes)


def merge_synteny_hashes(
    args,
    unambig_hashes=None,
    ambig_hashes=None,
    hasher=None,
    n_proteomes=None,
    file_writer=None,
):
    """Merge synteny hashes into proteomes."""
    hash_name = hasher.hash_name()
    plain_hash_name = hasher.hash_name(no_prefix=True)
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = syn.join(unambig_hashes, on=hash_name)
    syn = syn.join(ambig_hashes, on=hash_name)
    homology = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
    syn = pd.concat([homology, syn], axis=1)
    n_proteins = len(syn)
    syn["tmp.i"] = range(len(syn))
    shingled_vars = {
        plain_hash_name: np.array([np.nan] * n_proteins),
        "direction": np.array([""] * n_proteins),
        "self_count": np.array([np.nan] * n_proteins),
        "footprint": np.array([np.nan] * n_proteins),
        "frag_count": np.array([np.nan] * n_proteins),
        "ortho_count": np.array([np.nan] * n_proteins),
    }
    anchor_blocks = np.array([np.nan] * n_proteins)
    for hash_val, subframe in syn.groupby(by=["tmp.base"]):
        for unused_i, row in subframe.iterrows():
            footprint = row["tmp.footprint"]
            row_no = row["tmp.i"]
            anchor_vec = hasher.shingle(
                syn["hom.cluster"][row_no : row_no + footprint],
                row["tmp.base"],
                row["tmp.direction"],
                row[hash_name],
            )
            anchor_blocks[row_no : row_no + footprint] = anchor_vec
            for key in shingled_vars:
                shingled_vars[key][row_no : row_no + footprint] = row[
                    "tmp." + key
                ]
    syn["syn.anchor_id"] = pd.array(anchor_blocks, dtype=pd.UInt32Dtype())
    syn["syn.direction"] = pd.Categorical(
        shingled_vars["direction"], dtype=DIRECTIONAL_CATEGORY
    )
    syn["syn." + plain_hash_name] = pd.array(
        shingled_vars[plain_hash_name], dtype=pd.Int64Dtype()
    )
    del shingled_vars["direction"], shingled_vars[plain_hash_name]
    for key in shingled_vars:
        syn["syn." + key] = pd.array(
            shingled_vars[key], dtype=pd.UInt32Dtype()
        )
    write_tsv_or_parquet(
        syn, outpath / SYNTENY_FILE,
    )
    for cluster_id, subframe in syn.groupby(by=["hom.cluster"]):
        with file_writer(cluster_id) as fh:
            subframe[CLUSTER_COLS].dropna().to_csv(fh, header=False, sep="\t")
    in_synteny = n_proteins - syn["syn.anchor_id"].isnull().sum()
    n_assigned = n_proteins - syn["hom.cluster"].isnull().sum()
    # Do histogram of blocks
    anchor_counts = syn["syn.anchor_id"].value_counts()
    anchor_hist = pd.DataFrame(anchor_counts.value_counts()).sort_index()
    anchor_hist = anchor_hist.rename(columns={"syn.anchor_id": "self_count"})
    anchor_hist["pct_anchors"] = (
        anchor_hist["self_count"] * anchor_hist.index * 100.0 / n_assigned
    )
    write_tsv_or_parquet(anchor_hist, outpath / ANCHOR_HIST_FILE)
    ambig_anchors = 0
    ambig_proteins = 0
    for unused_anchor_id, subframe in syn.groupby(by=["syn.anchor_id"]):
        self_count = subframe["syn.self_count"].iloc[0]
        if self_count > 1:
            ambig_anchors += 1
            ambig_proteins += self_count
    # ambig = (syn["syn.self_count"] != 1).sum()
    avg_ortho = syn["syn.ortho_count"].mean()
    synteny_pct = in_synteny * 100.0 / n_assigned
    unambig_pct = (in_synteny - ambig_proteins) * 100.0 / n_assigned
    synteny_stats = {
        "idx": idx,
        "path": dotpath,
        "hom.clusters": n_assigned,
        "syn.anchors": in_synteny,
        "syn.ambig": ambig_anchors,
        "syn.assgn_pct": synteny_pct,
        "syn.unamb_pct": unambig_pct,
        "syn.fom": avg_ortho * 100.0 / n_proteomes,
    }
    return synteny_stats


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
