# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import sys
from itertools import combinations

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import dask.bag as db
import networkx as nx
import numpy as np
import pandas as pd
from dask.diagnostics import ProgressBar

# module imports
from .common import AMBIGUOUS_CODE
from .common import ANCHOR_HIST_FILE
from .common import ANCHORS_FILE
from .common import CLUSTERS_FILE
from .common import CLUSTERSYN_FILE
from .common import CODE_DICT
from .common import DISAMBIGUATED_CODE
from .common import HOMOLOGY_FILE
from .common import INDIRECT_CODE
from .common import LOCALLY_UNAMBIGUOUS_CODE
from .common import NON_AMBIGUOUS_CODE
from .common import PROTEOMOLOGY_FILE
from .common import PROTEOSYN_FILE
from .common import SPINNER_UPDATE_PERIOD
from .common import SYNTENY_FILE
from .common import SYNTENY_FILETYPE
from .common import UNAMBIGUOUS_CODE
from .common import calculate_adjacency_group
from .common import dotpath_to_path
from .common import hash_array
from .common import log_and_add_to_stats
from .common import logger
from .common import read_tsv_or_parquet
from .common import write_tsv_or_parquet
from .hash import SyntenyBlockHasher
from .mailboxes import DataMailboxes
from .mailboxes import ExternalMerge
from .merger import AmbiguousMerger

# global constants
__ALL__ = ["synteny_anchors"]
CLUSTER_COLS = ["syn.anchor.id", "syn.anchor.count", "syn.anchor.direction"]
JOIN_COLS = [
    "member_ids",
    "syn.anchor.sub_id",
    "syn.anchor.id",
    "syn.anchor.count",
    "syn.code",
    "frag.idx",
]
ANCHOR_COLS = [
    "path",
    "syn.code",
    "syn.anchor.count",
    "hom.cluster",
    "frag.id",
    "frag.pos",
    "hom.cl_size",
    "frag.direction",
    "frag.idx",
    "frag.is_chr",
    "frag.is_plas",
    "frag.is_scaf",
    "frag.prot_count",
    "frag.start",
    "prot.len",
    "prot.m_start",
    "prot.n_ambig",
    "prot.no_stop",
]
MAILBOX_SUBDIR = "mailboxes"

# CLI function
def synteny_anchors(
    k,
    peatmer,
    setname,
    click_loguru=None,
    write_ambiguous=True,
    thorny=True,
    disambig_adj_only=True,
):
    """Calculate synteny anchors."""
    #
    # Marshal input arguments
    #
    if k < 2:
        logger.error("k must be at least 2.")
        sys.exit(1)
    options = click_loguru.get_global_options()
    user_options = click_loguru.get_user_global_options()
    set_path = Path(setname)
    file_stats_path = set_path / PROTEOMOLOGY_FILE
    proteomes = read_tsv_or_parquet(file_stats_path)
    n_proteomes = len(proteomes)
    clusters = read_tsv_or_parquet(set_path / CLUSTERS_FILE)
    n_clusters = len(clusters)
    hasher = SyntenyBlockHasher(
        k=k,
        peatmer=peatmer,
        thorny=thorny,
        disambig_adj_only=disambig_adj_only,
    )
    logger.info(
        f"Calculating {hasher.hash_name(no_prefix=True)} synteny anchors"
        + f" for {n_proteomes} proteomes"
    )
    # durable argument list for passes
    arg_list = [
        (
            idx,
            row["path"],
        )
        for idx, row in proteomes.iterrows()
    ]
    runner = PassRunner(
        {
            "n_proteomes": n_proteomes,
            "set_path": set_path,
            "hasher": hasher,
            "quiet": options.quiet,
            "parallel": user_options["parallel"],
            "bag": db.from_sequence(arg_list),
            "merge_args": arg_list,
            "click_loguru": click_loguru,
        }
    )
    #
    # Make first three passes:
    #  1. hash and find unambiguous anchors
    #  2. disambiguate ambiguous anchors adjacent to unambiguous ones
    #  3. find non-ambiguous hashes uncovered by 2)
    #
    for pass_code in [
        UNAMBIGUOUS_CODE,
        DISAMBIGUATED_CODE,
        NON_AMBIGUOUS_CODE,
    ]:
        proteomes = runner.make_pass(pass_code, proteomes)
    runner.add_ambig_to_total_assigned()
    n_anchors = runner.get_total_assigned()
    #
    # Fourth pass -- merge, write anchor and homology info
    #
    join_mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / MAILBOX_SUBDIR / "join"),
        file_extension="tsv",
    )
    join_mb.write_tsv_headers(JOIN_COLS)
    cluster_mb = DataMailboxes(
        n_boxes=n_clusters,
        mb_dir_path=(set_path / MAILBOX_SUBDIR / "clusters"),
        file_extension="tsv",
    )
    cluster_mb.write_tsv_headers(CLUSTER_COLS)
    anchor_mb = DataMailboxes(
        n_boxes=n_anchors,
        mb_dir_path=(set_path / MAILBOX_SUBDIR / "anchors"),
        file_extension="tsv",
    )
    anchor_mb.write_tsv_headers(ANCHOR_COLS)
    proteomes = runner.make_pass(
        INDIRECT_CODE,
        proteomes,
        extra_kwargs={
            "join_mb": join_mb,
            "cluster_mb": cluster_mb,
            "anchor_mb": anchor_mb,
            "n_proteomes": n_proteomes,
            "write_ambiguous": write_ambiguous,
        },
    )
    write_tsv_or_parquet(
        proteomes, set_path / PROTEOSYN_FILE, remove_tmp=False
    )
    adjacency_stats = anchors_to_adjacency(
        set_path, n_proteomes, join_mb.open_then_delete
    )
    logger.info(f"adjacency_stats: {adjacency_stats}")
    return
    #
    # Write anchors
    #
    click_loguru.elapsed_time("Anchor writing")
    arg_list = [(i,) for i in range(n_anchors)]
    anchor_path = set_path / "synteny"
    anchor_path.mkdir(exist_ok=True)
    logger.info(f"Writing {n_anchors} synteny anchors to {anchor_path}:")
    if not options.quiet:
        ProgressBar(dt=SPINNER_UPDATE_PERIOD).register()
    if parallel:
        bag = db.from_sequence(arg_list)
        anchor_stats = bag.map(
            write_anchor,
            mailbox_reader=anchor_mb.open_then_delete,
            synteny_parent=anchor_path,
        ).compute()
    else:
        anchor_stats = []
        for args in arg_list:
            anchor_stats.append(
                write_anchor(
                    args,
                    mailbox_reader=anchor_mb.open_then_delete,
                    synteny_parent=anchor_path,
                )
            )
    anchor_mb.delete()
    anchor_stat_list = []
    for results in anchor_stats:
        if results is not None:
            anchor_stat_list += results
    anchor_frame = pd.DataFrame.from_dict(results)
    write_tsv_or_parquet(
        anchor_frame,
        set_path / ANCHORS_FILE,
        sort_cols=False,
    )
    #
    # Merge synteny into clusters
    #
    arg_list = [(i,) for i in range(n_clusters)]
    click_loguru.elapsed_time("Synteny joining")
    homology_path = set_path / "homology"
    logger.info(
        f"Joining synteny info to {n_clusters} clusters in {homology_path}:"
    )
    if not options.quiet:
        ProgressBar(dt=SPINNER_UPDATE_PERIOD).register()
    if parallel:
        bag = db.from_sequence(arg_list)
        cluster_stats = bag.map(
            join_synteny_to_clusters,
            mailbox_reader=cluster_mb.open_then_delete,
            cluster_parent=homology_path,
        ).compute()
    else:
        cluster_stats = []
        for args in arg_list:
            cluster_stats.append(
                join_synteny_to_clusters(
                    args,
                    mailbox_reader=cluster_mb.open_then_delete,
                    cluster_parent=homology_path,
                )
            )
    cluster_mb.delete()
    cluster_frame = pd.DataFrame.from_dict(cluster_stats)
    cluster_frame.set_index(["clust_id"], inplace=True)
    cluster_frame.sort_index(inplace=True)
    clusters = _concat_without_overlap(clusters, cluster_frame)
    write_tsv_or_parquet(
        clusters, set_path / CLUSTERSYN_FILE, float_format="%5.2f"
    )
    mean_gene_synteny = (
        clusters["in_synteny"].sum() * 100.0 / clusters["size"].sum()
    )
    mean_clust_synteny = clusters["synteny_pct"].mean()
    logger.info(
        f"Mean anchor coverage: {mean_gene_synteny: .1f}% (on proteins)"
    )
    logger.info(
        f"Mean cluster anchor coverage: {mean_clust_synteny:.1f}% (on clusters)"
    )
    click_loguru.elapsed_time(None)


class PassRunner:
    """Run a pass over all proteomes."""

    def __init__(self, std_kwargs):
        """Save initial pass info"""
        self.std_kwargs = std_kwargs
        self.last_code = None
        self.pass_name = "Hashing"
        self.merger_kw_dict = {
            UNAMBIGUOUS_CODE: {
                "count_key": "syn.anchor.count",
                "ordinal_key": "syn.anchor.id",
                "ambig_ordinal_key": "tmp.ambig.id",
            },
            DISAMBIGUATED_CODE: {
                "count_key": "tmp.disambig.anchor.count",
                "ordinal_key": "tmp.disambig.anchor.id",
                "alt_hash": True,
            },
            NON_AMBIGUOUS_CODE: {
                "count_key": "tmp.nonambig.anchor.count",
                "ordinal_key": "tmp.nonambig.anchor.id",
                "ambig_count_key": "tmp.ambig.anchor.count",
                "ambig_ordinal_key": "tmp.ambig.anchor.id",
            },
        }
        self.merge_function_dict = {
            UNAMBIGUOUS_CODE: calculate_synteny_hashes,
            DISAMBIGUATED_CODE: merge_unambig_hashes,
            NON_AMBIGUOUS_CODE: merge_disambig_hashes,
            INDIRECT_CODE: merge_nonambig_hashes,
        }
        self.n_assigned_list = []
        self.ambig = None
        self.unambig = None
        self.log_ambig = False

    def make_pass(
        self,
        code,
        proteomes,
        extra_kwargs=None,
    ):
        """Make a calculate-merge pass over each proteome."""
        self.std_kwargs["click_loguru"].elapsed_time(self.pass_name)
        if self.unambig is not None:
            if self.log_ambig:
                ambig_msg = f" and {len(self.ambig)} ambiguous"
            else:
                ambig_msg = ""
            logger.info(
                f"Merging {len(self.unambig)} {CODE_DICT[self.last_code]}"
                + f"({self.last_code}){ambig_msg} synteny anchors into proteomes"
            )
        if extra_kwargs is None:
            extra_kwargs = {}
        kwargs = {
            "unambig": self.unambig,
            "ambig": self.ambig,
            "hasher": self.std_kwargs["hasher"],
        }
        if code in self.merger_kw_dict:
            mailboxes = DataMailboxes(
                n_boxes=self.std_kwargs["n_proteomes"],
                mb_dir_path=(
                    self.std_kwargs["set_path"]
                    / MAILBOX_SUBDIR
                    / CODE_DICT[code]
                ),
            )
            mailboxes.write_headers("hash\n")
            kwargs["mailboxes"] = mailboxes
        merge_func = self.merge_function_dict[code]
        if not self.std_kwargs["quiet"]:
            ProgressBar(dt=SPINNER_UPDATE_PERIOD).register()
        if self.std_kwargs["parallel"]:
            stats_list = (
                self.std_kwargs["bag"]
                .map(merge_func, **kwargs, **extra_kwargs)
                .compute()
            )
        else:
            stats_list = []
            for args in self.std_kwargs["merge_args"]:
                stats_list.append(merge_func(args, **kwargs, **extra_kwargs))
        stats = (
            pd.DataFrame.from_dict(stats_list).set_index("idx").sort_index()
        )
        proteomes = log_and_add_to_stats(proteomes, stats)
        if code in self.merger_kw_dict:
            merger = ExternalMerge(
                file_path_func=mailboxes.path_to_mailbox,
                n_merge=self.std_kwargs["n_proteomes"],
            )
            merger.init("hash")
            merge_counter = AmbiguousMerger(
                start_base=sum(self.n_assigned_list),
                **self.merger_kw_dict[code],
            )
            self.unambig, self.ambig = merger.merge(merge_counter)
            mailboxes.delete()
        self.n_assigned_list.append(len(self.unambig))
        self.last_code = code
        self.pass_name = CODE_DICT[code]
        return proteomes

    def add_ambig_to_total_assigned(self):
        """Include the most recent ambiguous assignment in the total."""
        self.n_assigned_list.append(len(self.ambig))
        self.log_ambig = True

    def get_total_assigned(self):
        """Return the total assigned anchors."""
        return sum(self.n_assigned_list)


def calculate_synteny_hashes(
    args, mailboxes=None, hasher=None, unambig=None, ambig=None
):
    """Calculate synteny hashes for proteins per-genome."""
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    hom = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
    hom["tmp.nan_group"] = (
        (hom["hom.cluster"].isnull()).astype(int).cumsum() + 1
    ) * (~hom["hom.cluster"].isnull())
    hom.replace(to_replace={"tmp.nan_group": 0}, value=pd.NA, inplace=True)
    hash_name = hasher.hash_name()
    syn_list = []
    if hasher.thorny:  # drop rows
        hom = hom[hom["hom.cluster"].notna()]
    for unused_id_tuple, subframe in hom.groupby(
        by=["frag.id", "tmp.nan_group"]
    ):
        syn_list.append(hasher.calculate(subframe["hom.cluster"]))
    del hom["tmp.nan_group"]
    syn = hom.join(
        pd.concat([df for df in syn_list if df is not None], axis=0)
    )
    del syn_list
    write_tsv_or_parquet(syn, outpath / SYNTENY_FILE, remove_tmp=False)
    syn["tmp.self_count"] = pd.array(
        syn[hash_name].map(syn[hash_name].value_counts()),
        dtype=pd.UInt32Dtype(),
    )
    unique_hashes = (
        syn[[hash_name, "tmp.self_count"]]
        .drop_duplicates(subset=[hash_name])
        .dropna(how="any")
    )
    unique_hashes = unique_hashes.set_index(hash_name).sort_index()
    with mailboxes.locked_open_for_write(idx) as file_handle:
        unique_hashes.to_csv(file_handle, header=False, sep="\t")
    return {
        "idx": idx,
        "path": dotpath,
        "hom.clusters": syn["hom.cluster"].notna().sum(),
        "syn.hashes.n": syn[hash_name].notna().sum(),
    }


def merge_unambig_hashes(
    args,
    unambig=None,
    ambig=None,
    hasher=None,
    mailboxes=None,
):
    """Merge unambiguous synteny hashes into proteomes per-proteome."""
    hash_name = hasher.hash_name()
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = _join_on_col_with_na(syn, unambig, hash_name)
    syn = _join_on_col_with_na(syn, ambig, hash_name)
    syn["syn.code"] = pd.NA
    syn["syn.code"] = _fill_col1_val_where_col2_notna(
        syn["syn.code"], syn["syn.anchor.id"], UNAMBIGUOUS_CODE
    )
    # Calculate disambiguation hashes and write them out for merge
    disambig_frame_list = []
    for unused_frag, subframe in syn.groupby(by=["frag.id"]):
        disambig_frame_list.append(hasher.calculate_disambig_hashes(subframe))
    disambig_fr = pd.concat(
        [df for df in disambig_frame_list if df is not None]
    )
    disambig_fr = disambig_fr.dropna(how="all")
    syn = syn.join(disambig_fr)
    write_tsv_or_parquet(syn, outpath / SYNTENY_FILE, remove_tmp=False)
    # Write out unified upstream/downstream hash values
    merged_hashes = pd.concat(
        [
            _rename_and_fill_alt(syn, "tmp.disambig.up", "tmp.disambig.down"),
            _rename_and_fill_alt(syn, "tmp.disambig.down", "tmp.disambig.up"),
        ],
        ignore_index=True,
    )
    merged_hashes["self_count"] = pd.array(
        merged_hashes["hash"].map(merged_hashes["hash"].value_counts()),
        dtype=pd.UInt32Dtype(),
    )
    merged_hashes = merged_hashes.reindex(
        columns=["hash", "self_count", "alt_hash"]
    )
    unique_hashes = (
        merged_hashes.drop_duplicates(subset=["hash"])
        .set_index("hash")
        .sort_index()
    )
    del merged_hashes
    with mailboxes.locked_open_for_write(idx) as file_handle:
        unique_hashes.to_csv(file_handle, header=False, sep="\t")
    return {
        "idx": idx,
        "path": dotpath,
        "syn.anchors.unambiguous": _count_code(
            syn["syn.code"], UNAMBIGUOUS_CODE
        ),
    }


def merge_disambig_hashes(
    args,
    unambig=None,
    ambig=None,
    hasher=None,
    mailboxes=None,
):
    """Merge disambiguated synteny hashes into proteomes per-proteome."""
    idx, dotpath = args
    plain_hash_name = hasher.hash_name(no_prefix=True)
    hash_name = "syn." + plain_hash_name
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = _join_on_col_with_na(syn, unambig, "tmp.disambig.up")
    syn = _join_on_col_with_na(syn, unambig, "tmp.disambig.down")
    for dup_col in [
        "tmp.disambig.anchor.count",
        "tmp.disambig.anchor.id",
    ]:
        xcol = dup_col + "_x"
        ycol = dup_col + "_y"
        syn[dup_col] = syn[xcol].fillna(syn[ycol])
        del syn[xcol], syn[ycol]
    syn["syn.anchor.id"] = syn["syn.anchor.id"].fillna(
        syn["tmp.disambig.anchor.id"]
    )
    syn["syn.anchor.count"] = syn["syn.anchor.count"].fillna(
        syn["tmp.disambig.anchor.count"]
    )
    syn["syn.code"] = _fill_col1_val_where_col2_notna(
        syn["syn.code"], syn["tmp.disambig.anchor.id"], DISAMBIGUATED_CODE
    )
    # Delete some non-needed tmp columns
    non_needed_cols = [
        "tmp.disambig.anchor.count",
        "tmp.disambig.anchor.id",
        "tmp.disambig.up",
        "tmp.disambig.down",
    ]
    syn = syn.drop(columns=non_needed_cols)
    # null hashes are already assigned
    syn[hash_name][syn["syn.anchor.id"].notna()] = pd.NA
    write_tsv_or_parquet(syn, outpath / SYNTENY_FILE, remove_tmp=False)
    # Write out non-ambiguous hashes
    syn["tmp.self_count"] = pd.array(
        syn[hash_name].map(syn[hash_name].value_counts()),
        dtype=pd.UInt32Dtype(),
    )
    unique_hashes = (
        syn[[hash_name, "tmp.self_count"]]
        .drop_duplicates(subset=[hash_name])
        .dropna(how="any")
    )
    unique_hashes = unique_hashes.set_index(hash_name).sort_index()
    with mailboxes.locked_open_for_write(idx) as file_handle:
        unique_hashes.to_csv(file_handle, header=False, sep="\t")
    # logger.debug(f"{dotpath} has {syn['syn.anchor.id'].notna().sum()} assignments")
    return {
        "idx": idx,
        "path": dotpath,
        "syn.anchors.disambiguated": _count_code(
            syn["syn.code"], DISAMBIGUATED_CODE
        ),
    }


def merge_nonambig_hashes(
    args,
    join_mb=None,
    unambig=None,
    ambig=None,
    hasher=None,
    n_proteomes=None,
    cluster_mb=None,
    anchor_mb=None,
    write_ambiguous=True,
):
    """Merge disambiguated synteny hashes into proteomes per-proteome."""
    idx, dotpath = args
    plain_hash_name = hasher.hash_name(no_prefix=True)
    hash_name = "syn." + plain_hash_name
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = _join_on_col_with_na(syn, unambig, hash_name)
    #
    # Do the indirects (formerly ambig made nonambig)
    #
    syn["syn.anchor.id"] = syn["syn.anchor.id"].fillna(
        syn["tmp.nonambig.anchor.id"]
    )
    syn["syn.anchor.count"] = syn["syn.anchor.count"].fillna(
        syn["tmp.nonambig.anchor.count"]
    )
    syn["syn.code"] = _fill_col1_val_where_col2_notna(
        syn["syn.code"], syn["tmp.nonambig.anchor.id"], INDIRECT_CODE
    )
    #
    # Do the nonambig (w.r.t. this proteome) and ambig, if requested
    #
    syn = _join_on_col_with_na(syn, ambig, hash_name)
    n_proteins = len(syn)
    syn["tmp.i"] = range(n_proteins)
    ambig_code = syn["syn.code"].copy()
    ambig_ids = syn["tmp.ambig.anchor.id"].copy()
    ambig_counts = syn["tmp.ambig.anchor.count"].copy()
    for unused_ambig_id, subframe in syn.groupby(by=["tmp.ambig.anchor.id"]):
        ambig_n = len(subframe)
        for unused_i, row in subframe.iterrows():
            row_no = row["tmp.i"]
            if ambig_n == 1:
                ambig_code.iloc[row_no] = LOCALLY_UNAMBIGUOUS_CODE
            elif ambig_n > 1:
                if write_ambiguous:
                    ambig_code.iloc[row_no] = AMBIGUOUS_CODE
                else:
                    ambig_ids.iloc[row_no] = pd.NA
                    ambig_counts.iloc[row_no] = pd.NA
    syn["syn.anchor.id"] = syn["syn.anchor.id"].fillna(ambig_ids)
    syn["syn.anchor.count"] = syn["syn.anchor.count"].fillna(ambig_counts)
    syn["syn.code"] = syn["syn.code"].fillna(ambig_code)
    del ambig_code, ambig_ids, ambig_counts
    #
    # Hsh footprint and direction are anchor properties, where set
    #
    syn = syn.rename(
        columns={
            "syn.hash.footprint": "syn.anchor.footprint",
            "syn.hash.direction": "syn.anchor.direction",
        }
    )
    n_anchors = syn["syn.anchor.id"].notna().sum()  # before shingling
    #
    # Do shingling
    #
    # shingle_id = np.array([np.nan] * n_proteins)
    # shingle_count = np.array([np.nan] * n_proteins)
    # shingle_code = syn["syn.code"].to_numpy()
    # shingle_direction = syn["syn.anchor.direction"].to_numpy()
    # shingle_sub = np.array([np.nan] * n_proteins)
    with join_mb.locked_open_for_write(idx) as file_handle:
        for anchor_count, subframe in syn.groupby(by=["syn.anchor.count"]):
            for unused_i, row in subframe.iterrows():
                anchor_id = row["syn.anchor.id"]
                if pd.isnull(anchor_id):
                    continue
                first_row = row["tmp.i"]
                last_row = first_row + row["syn.anchor.footprint"]
                # shingle_id[first_row:last_row] = anchor_id
                # shingle_count[first_row:last_row] = anchor_count
                # shingle_code[first_row:last_row] = row["syn.code"]
                # shingle_direction[first_row:last_row] = row[
                #    "syn.anchor.direction"
                # ]
                # shingle_sub[first_row:last_row] = hasher.shingle(
                #    syn["hom.cluster"][first_row:last_row],
                #    row["syn.anchor.direction"],
                #    row[hash_name],
                # )
                shingle_fr = pd.DataFrame(
                    {
                        "member_ids": syn.iloc[first_row:last_row].index,
                        "syn.anchor.sub_id": hasher.shingle(
                            syn["hom.cluster"][first_row:last_row],
                            row["syn.anchor.direction"],
                            row[hash_name],
                        ),
                    }
                )
                shingle_fr["syn.anchor.id"] = anchor_id
                shingle_fr["syn.anchor.count"] = anchor_count
                shingle_fr["syn.code"] = row["syn.code"]
                shingle_fr["frag.idx"] = row["frag.idx"]
                shingle_fr.to_csv(file_handle, header=False, sep="\t")
    # syn["syn.anchor.id"] = shingle_id
    # syn["syn.anchor.count"] = shingle_count
    # syn["syn.code"] = shingle_code
    # syn["syn.anchor.direction"] = shingle_direction
    # syn["syn.anchor.sub_id"] = shingle_sub
    # del shingle_id, shingle_count, shingle_code, shingle_sub
    # Delete non-needed (but non-tmp) columns
    # non_needed_cols = [hash_name]
    # syn = syn.drop(columns=non_needed_cols)
    # Write proteome file
    # write_tsv_or_parquet(
    #    syn,
    #    outpath / SYNTENY_FILE,
    # )
    # Write anchor info to mailbox
    # for anchor_id, subframe in syn.groupby(by=["syn.anchor.id"]):
    #    anchor_frame = subframe.copy()
    #    anchor_frame["path"] = dotpath
    #    with anchor_mb.locked_open_for_write(anchor_id) as file_handle:
    #        anchor_frame[ANCHOR_COLS].to_csv(
    #            file_handle, header=False, sep="\t"
    #        )
    # for cluster_id, subframe in syn.groupby(by=["hom.cluster"]):
    #    with cluster_mb.locked_open_for_write(cluster_id) as file_handle:
    #        subframe[CLUSTER_COLS].to_csv(file_handle, header=False, sep="\t")
    in_synteny = syn["syn.anchor.id"].notna().sum()
    n_assigned = syn["hom.cluster"].notna().sum()
    avg_ortho = syn["syn.anchor.count"].mean()
    synteny_pct = in_synteny * 100.0 / n_assigned
    n_ambig = _count_code(syn["syn.code"], AMBIGUOUS_CODE)
    n_nonambig = in_synteny - n_ambig
    nonambig_pct = n_nonambig * 100.0 / n_assigned
    synteny_stats = {
        "idx": idx,
        "path": dotpath,
        "syn.anchors.indirect_unambiguous": _count_code(
            syn["syn.code"], INDIRECT_CODE
        ),
        "syn.anchors.locally_unambiguous": _count_code(
            syn["syn.code"], LOCALLY_UNAMBIGUOUS_CODE
        ),
        "syn.anchors.ambiguous": n_ambig,
        "syn.anchors.nonambiguous": n_nonambig,
        "syn.anchors.nonambig_pct": nonambig_pct,
        "syn.anchors.base": n_anchors,
        "syn.anchors.total": in_synteny,
        "syn.anchors.total_pct": synteny_pct,
        "syn.orthogenomic_pct": avg_ortho * 100.0 / n_proteomes,
    }
    return synteny_stats


def join_synteny_to_clusters(args, cluster_parent=None, mailbox_reader=None):
    """Read homology info from mailbox and join it to proteome file."""
    idx = args[0]
    cluster_path = cluster_parent / f"{idx}.parq"
    cluster = pd.read_parquet(cluster_path)
    n_cluster = len(cluster)
    with mailbox_reader(idx) as file_handle:
        synteny_frame = pd.read_csv(
            file_handle, sep="\t", index_col=0
        ).convert_dtypes()
        in_synteny = len(synteny_frame)
    # cluster files are unusual in that I don't bother to version them,
    # so overlapping info has to be deleted each time
    clust_syn = _concat_without_overlap(cluster, synteny_frame)
    write_tsv_or_parquet(clust_syn, cluster_path)
    anchor_count = clust_syn["syn.anchor.id"].value_counts()
    anchor_frag_counts = [0]
    for unused_id_tuple, subframe in clust_syn.groupby(
        by=["syn.anchor.id", "path"]
    ):
        if len(subframe) == 1:
            anchor_frag_counts.append(1)
        else:
            anchor_frag_counts.append(len(subframe["frag.idx"].value_counts()))
    return {
        "clust_id": idx,
        "in_synteny": in_synteny,
        "n_anchors": len(anchor_count),
        "max_frags_per_anch": max(anchor_frag_counts),
        "synteny_pct": in_synteny * 100.0 / n_cluster,
    }


def write_anchor(args, synteny_parent=None, mailbox_reader=None):
    """Read synteny anchor info from mailbox and join it to synteny file."""
    idx = args[0]
    with mailbox_reader(idx) as file_handle:
        anchor_frame = pd.read_csv(
            file_handle, sep="\t", index_col=0
        ).convert_dtypes()
    in_anchor = len(anchor_frame)
    if in_anchor == 0:
        return None
    # drop any duplicated ID's--normally shouldn't happen
    anchor_frame.drop(
        anchor_frame[anchor_frame.index.duplicated()].index, inplace=True
    )
    anchor_frame.sort_values(
        by=["syn.anchor.sub_id", "frag.idx", "frag.pos"], inplace=True
    )
    # Make a dictionary of common anchor properties, order will be kept
    anchor_props = {
        "anchor.id": idx,
        "sub": None,
        "code": None,
        "count": None,
        "n": None,
        "n_ambig": None,
        "n_adj": None,
        "adj_groups": None,
        "frag.direction": None,
        "syn.anchor.direction": None,
        "anchor.subframe.ok": True,
        "hash": None,
    }
    code_set = set(anchor_frame["syn.code"])
    for test_code in CODE_DICT.keys():
        if test_code in code_set:
            anchor_props["code"] = test_code
            break
    bad_subframe = False
    prop_list = []
    for sub_no, subframe in anchor_frame.groupby(by=["syn.anchor.sub_id"]):
        (subanchor_props, anchor_subframe, bad_subframe) = _subframe_props(
            anchor_props, subframe, sub_no
        )
        if bad_subframe:
            break
        write_tsv_or_parquet(
            anchor_subframe,
            synteny_parent / f"{idx}.{sub_no}.{SYNTENY_FILETYPE}",
            sort_cols=False,
        )
        prop_list.append(subanchor_props)
    if bad_subframe:  # Probably means a hash collision
        logger.error(f"bad anchor set {idx}")
        prop_list = []
        sub_no = 0
        anchor_props["anchor.subframe.ok"] = False
        for cluster_id, subframe in anchor_frame.groupby(by=["hom.cluster"]):
            (
                subanchor_props,
                anchor_subframe,
                unused_bad_subframe,
            ) = _subframe_props(anchor_props, subframe, sub_no)
            write_tsv_or_parquet(
                anchor_subframe,
                synteny_parent / f"{idx}.{sub_no}.{SYNTENY_FILETYPE}",
                sort_cols=False,
            )
            sub_no += 1
            prop_list.append(subanchor_props)
    return prop_list


def _subframe_props(anchor_props, subframe, sub_no):
    """Calculate subframe properties and write subframe"""
    anchor_subframe = subframe.copy()
    subanchor_props = anchor_props.copy()
    subanchor_props["sub"] = sub_no
    anchor_dir_set = set(anchor_subframe["syn.anchor.direction"])
    if len(anchor_dir_set) == 1:
        subanchor_props["syn.anchor.direction"] = list(anchor_dir_set)[0]
    frag_dir_set = set(anchor_subframe["frag.direction"])
    if len(frag_dir_set) == 1:
        subanchor_props["frag.direction"] = list(frag_dir_set)[0]
    subanchor_props["count"] = anchor_subframe["syn.anchor.count"].iloc[0]
    subanchor_props["n_ambig"] = _count_code(
        anchor_subframe["syn.code"], AMBIGUOUS_CODE
    )
    hom_clust_set = set(anchor_subframe["hom.cluster"])
    if len(hom_clust_set) == 1:
        subanchor_props[f"anchor.{sub_no}.cluster"] = list(hom_clust_set)[0]
    else:
        bad_subframe = True
    del (
        anchor_subframe["syn.anchor.count"],
        anchor_subframe["syn.anchor.sub_id"],
    )
    subanchor_props["n"] = len(anchor_subframe)
    subanchor_props["hash"] = hash_array(
        np.sort(anchor_subframe.index.to_numpy())
    )
    (
        subanchor_props["n_adj"],
        subanchor_props["adj_groups"],
        unused_adj_group,
    ) = calculate_adjacency_group(
        anchor_subframe["frag.pos"], anchor_subframe["frag.idx"]
    )
    return subanchor_props, anchor_subframe, bad_subframe


def _concat_without_overlap(df1, df2):
    """Concatenate df2 on top of df1."""
    overlapping = set(df1.columns).intersection(df2.columns)
    if len(overlapping) > 0:
        df1 = df1.drop(columns=overlapping)
    return pd.concat([df1, df2], axis=1)


def _rename_and_fill_alt(df1, key, alt_key):
    """Rename columns and zero-fill alternate."""
    df2 = df1[[key]].rename(columns={key: "hash"})
    df2["alt_hash"] = df1[alt_key].fillna(0)
    return df2.dropna(how="any")


def _join_on_col_with_na(left, right, col_name):
    """Join on a temporary column of type 'O'."""
    tmp_col_name = "tmp." + col_name
    left[tmp_col_name] = left[col_name].astype("O")
    merged = pd.merge(
        left, right, left_on=tmp_col_name, right_index=True, how="left"
    )
    del merged[tmp_col_name]
    return merged


def _fill_col1_val_where_col2_notna(col1, col2, val):
    """Set col1 to val  where col2 is not NA if col1 is not set."""
    fill_ser = col1.copy()
    fill_ser[col2.notna()] = val
    return col1.fillna(fill_ser)


def _count_code(code_ser, code):
    """Counts number of occurrances of code in code_ser."""
    return (code_ser == code).sum()


def anchors_to_adjacency(set_path, n_proteomes, mailbox_reader):
    """Merge adjacencies and produce and adjacency graph."""
    frame_list = []
    for idx in range(n_proteomes):
        with mailbox_reader(idx) as file_handle:
            frame_list.append(
                pd.read_csv(
                    file_handle, sep="\t", index_col=0
                ).convert_dtypes()
            )
    nodes = pd.concat(
        frame_list,
        ignore_index=True,
    )
    del frame_list
    graph = nx.Graph()
    for unused_tuple, subframe in nodes.groupby(
        by=["syn.anchor.id", "syn.anchor.sub_id"]
    ):
        ids = subframe["member_ids"]
        n_ids = len(ids)
        graph.add_nodes_from(ids)
        if n_ids > 1:
            edges = combinations(ids, 2)
            graph.add_edges_from(edges, weight=n_ids)
    outpath = set_path / ANCHORS_FILE
    summarypath = outpath.parent / (
        outpath.name[: -len(outpath.suffix)] + "_summary.tsv"
    )
    histpath = outpath.parent / (
        outpath.name[: -len(outpath.suffix)] + "_hist.tsv"
    )
    components = [
        c
        for c in sorted(nx.connected_components(graph), key=len, reverse=True)
        if len(c) > 1
    ]
    fh = outpath.open("w")
    fh.write("idx\tcluster_id\tsize\tmembers\n")
    n_items = 0
    count_list = []
    hash_list = []
    id_list = []
    for i, comp in enumerate(components):
        component = np.sort(pd.Index(list(comp)).to_numpy())
        id_list.append(i)
        size = len(comp)
        count_list.append(size)
        hash_list.append(hash_array(component))
        for node in component:
            fh.write(f"{n_items}\t{i}\t{size}\t{node}\n")
            n_items += 1
    fh.close()
    n_clusts = len(count_list)
    del graph, components
    cluster_counts = pd.DataFrame({"size": count_list})
    largest_cluster = cluster_counts["size"].max()
    cluster_hist = (
        pd.DataFrame(cluster_counts.value_counts()).sort_index().reset_index()
    )
    cluster_hist = cluster_hist.set_index("size")
    cluster_hist = cluster_hist.rename(columns={0: "n"})
    cluster_hist["item_pct"] = (
        cluster_hist["n"] * cluster_hist.index * 100.0 / n_items
    )
    cluster_hist.to_csv(histpath, sep="\t", float_format="%5.2f")
    cluster_hist["cluster_pct"] = cluster_hist["n"] * 100.0 / n_clusts
    cluster_hist.to_csv(histpath, sep="\t", float_format="%5.2f")
    clusters = pd.DataFrame(
        {"anchor.id": id_list, "count": count_list, "hash": hash_list}
    )
    clusters.to_csv(summarypath, sep="\t")
    stats_dict = {
        "in_anchor": n_items,
        "syn.anchors.n": n_clusts,
        "syn.anchors.largest": largest_cluster,
    }
    return stats_dict


def intersect_anchors(set1_file, set2_file):
    set1_path = Path(set1_file)
    set2_path = Path(set2_file)
    set1_fr = pd.read_csv(set1_path, sep="\t", index_col=0)
    set2_fr = pd.read_csv(set2_path, sep="\t", index_col=0)
    set1_dict = {}
    for cluster_id, subframe in set1_fr.groupby(by=["cluster_id"]):
        set1_dict[cluster_id] = set(subframe["members"].to_numpy())
    set2_dict = {}
    for cluster_id, subframe in set2_fr.groupby(by=["cluster_id"]):
        set2_dict[cluster_id] = set(subframe["members"].to_numpy())
    identity_sets = []
    s1_subset = []
    s2_subset = []
    incongruent = []
    match_keys = list(set2_dict.keys())
    for key1 in set1_dict:
        s1 = set1_dict[key1]
        for i, key2 in enumerate(match_keys):
            s2 = set2_dict[key2]
            if len(s1.intersection(s2)) == 0:
                continue
            elif s1 == s2:
                identity_sets.append(
                    (
                        key1,
                        key2,
                    )
                )
                match_keys.pop(i)
                break
            elif s1.issubset(s2):
                s1_subset.append(
                    (
                        key1,
                        key2,
                    )
                )
                match_keys.pop(i)
                break
            elif s2.issubset(s1):
                s2_subset.append(
                    (
                        key1,
                        key2,
                    )
                )
                match_keys.pop(i)
                break
            else:
                incongruent.append(
                    (
                        key1,
                        key2,
                    )
                )
    logger.info(f"set 1 ({set1_file}): {len(set1_dict)}")
    logger.info(f"set 2 ({set2_file}): {len(set2_dict)}")
    min_sets = min(len(set1_dict), len(set2_dict))
    id_len = len(identity_sets)
    id_pct = id_len * 100.0 / min_sets
    logger.info(f"identity: {id_len} ({id_pct:.1f}%)")
    s1_len = len(s1_subset)
    s1_pct = s1_len * 100.0 / min_sets
    logger.info(f"set 1 is subset: {s1_len} ({s1_pct:.1f}%)")
    s2_len = len(s2_subset)
    s2_pct = s2_len * 100.0 / min_sets
    logger.info(f"set 2 is subset: {s2_len} ({s2_pct:.1f}%)")
    incon_len = len(incongruent)
    incon_pct = incon_len * 100.0 / min_sets
    logger.info(f"incongruent: {incon_len}({incon_pct}%)")
