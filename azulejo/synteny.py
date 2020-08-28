# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import array
import sys

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import click
import dask.bag as db
import numpy as np
import pandas as pd
import sh
from dask.diagnostics import ProgressBar
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import ANCHOR_HIST_FILE
from .common import CLUSTERS_FILE
from .common import CLUSTERSYN_FILE
from .common import HOMOLOGY_FILE
from .common import PROTEOMOLOGY_FILE
from .common import PROTEOSYN_FILE
from .common import SYNTENY_FILE
from .common import dotpath_to_path
from .common import log_and_add_to_stats
from .common import read_tsv_or_parquet
from .common import write_tsv_or_parquet
from .mailboxes import DataMailboxes
from .mailboxes import ExternalMerge
from .merger import AmbiguousMerger
from .hash import SyntenyBlockHasher

# global constants
CLUSTER_COLS = [
    "syn.anchor_id",
]
MERGE_COLS = ["syn.hash.self_count", "frag.idx"]
REMERGE_COLS = ["tmp.disambig.self_count", "frag.idx"]
DEFAULT_K = 2


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
    "--nonambig/--no-nonambig",
    is_flag=True,
    default=True,
    show_default=True,
    help="Fill with non-ambiguous hashes.",
)
@click.option(
    "--shingle/--no-shingle",
    is_flag=True,
    default=True,
    show_default=True,
    help="Shingle non-ambiguous hashes over k-mer.",
)
@click.argument("setname")
def synteny_anchors(k, peatmer, setname, greedy, nonambig, shingle):
    """Calculate synteny anchors."""
    if k < 2:
        logger.error("k must be at least 2.")
        sys.exit(1)
    options = click_loguru.get_global_options()
    user_options = click_loguru.get_user_global_options()
    parallel = user_options["parallel"]
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
    logger.info(
        f"Unambiguous anchor IDs in range {unambig_hashes['tmp.base'].min()}"
        + f" to {unambig_hashes['tmp.base'].max()+k}"
    )
    logger.info(
        f"Ambiguous anchor IDs in range {ambig_hashes['tmp.base.ambig'].min()}"
        + f" to {ambig_hashes['tmp.base.ambig'].max()+k}"
    )
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
            shingle=shingle,
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
                    shingle=shingle,
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
        f"Reducing {proteomes['tmp.hashes.pending_disambig'].sum()} disambiguation hashes"
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
        graph_path=set_path / "synteny.gml",
        start_base=last_anchor_base,
        k=k,
        graph=fragment_synteny_graph,
    )
    (
        disambig_hashes,
        still_ambig_hashes,
        fragment_synteny_graph,
    ) = disambig_merger.merge(disambig_merge_counter)

    logger.info(
        f"Disambiguated anchor IDs in range {disambig_hashes['tmp.disambig.base'].min()}"
        + f" to {disambig_hashes['tmp.disambig.base'].max()+k}"
    )
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
            + f"and {len(still_ambig_hashes)} still-ambiguous synteny anchors into {n_proteomes}"
            " proteomes"
        )
        if nonambig:
            logger.info(
                "Non-ambiguous instances of ambiguous hashes will be used to fill."
            )
        ProgressBar().register()
    if parallel:
        disambig_stats_list = bag.map(
            merge_disambig_hashes,
            disambig_hashes=disambig_hashes,
            still_ambig_hashes=still_ambig_hashes,
            hasher=hasher,
            n_proteomes=n_proteomes,
            cluster_writer=cluster_mb.locked_open_for_write,
            greedy=greedy,
            nonambig=nonambig,
            shingle=shingle,
        ).compute()
    else:
        for args in arg_list:
            disambig_stats_list.append(
                merge_disambig_hashes(
                    args,
                    disambig_hashes=disambig_hashes,
                    still_ambig_hashes=still_ambig_hashes,
                    hasher=hasher,
                    n_proteomes=n_proteomes,
                    cluster_writer=cluster_mb.locked_open_for_write,
                    greedy=greedy,
                    nonambig=nonambig,
                    shingle=shingle,
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
    proteomes = concat_without_overlap(clusters, anchor_frame)
    write_tsv_or_parquet(
        proteomes, set_path / CLUSTERSYN_FILE, float_format="%5.2f"
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
    """Concatenate df2 on top of df1."""
    overlapping = set(df1.columns).intersection(df2.columns)
    if len(overlapping) > 0:
        df1 = df1.drop(columns=overlapping)
    return pd.concat([df1, df2], axis=1)


def calculate_synteny_hashes(args, mailboxes=None, hasher=None):
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
    for unused_id_tuple, subframe in hom.groupby(
        by=["frag.id", "tmp.nan_group"]
    ):
        syn_list.append(hasher.calculate(subframe["hom.cluster"]))
    syn = hom.join(
        pd.concat([df for df in syn_list if df is not None], axis=0)
    )
    del syn_list
    syn["tmp.i"] = pd.array(range(len(syn)), dtype=pd.UInt32Dtype())
    syn["syn.hash.self_count"] = pd.array(
        syn[hash_name].map(syn[hash_name].value_counts()),
        dtype=pd.UInt32Dtype(),
    )
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
        "syn.hashes.n": in_hash,
        "syn.hash_pct": hash_pct,
    }
    return hash_stats


def merge_unambig_hashes(
    args,
    unambig_hashes=None,
    ambig_hashes=None,
    hasher=None,
    ambig_mb=None,
    shingle=True,
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
    if shingle:
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
    else:
        syn["syn.anchor_id"] = syn["tmp.base"]
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
    syn["tmp.disambig.self_count"] = pd.array(
        syn["syn.disambig_upstr"].map(
            syn["syn.disambig_upstr"].value_counts()
        ),
        dtype=pd.UInt32Dtype(),
    )
    upstream_hashes = (
        syn[["syn.disambig_upstr"] + REMERGE_COLS]
        .dropna(how="any")
        .rename(columns={"syn.disambig_upstr": "hash"})
    )
    syn["tmp.disambig.self_count"] = pd.array(
        syn["syn.disambig_downstr"].map(
            syn["syn.disambig_downstr"].value_counts()
        ),
        dtype=pd.UInt32Dtype(),
    )
    downstream_hashes = (
        syn[["syn.disambig_downstr"] + REMERGE_COLS]
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
    hasher=None,
    n_proteomes=None,
    cluster_writer=None,
    greedy=False,
    nonambig=False,
    shingle=True,
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
    if shingle:
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
        disambig_ser = pd.Series(disambig_fills, index=syn.index)
        syn["syn.anchor_id"] = syn["syn.anchor_id"].fillna(
            pd.Series(disambig_fills, index=syn.index)
        )
    else:
        disambig_ser = syn["tmp.disambig.base"]
    syn["syn.anchor_id"] = syn["syn.anchor_id"].fillna(disambig_ser)
    n_disambiguated = (syn["syn.anchor_id"] == disambig_ser).sum()
    # Deal with ambiguous hashes by adding non-ambiguous examples
    if nonambig:
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
    else:
        n_nonambig = 0
    # Delete temporary columns
    non_needed_cols = [
        "tmp.i",
        "syn.hash.disambig.ortho_count.ambig",
        "syn.hash.ortho_count.ambig",
        "syn.hash.footprint",
        # "syn.hash.direction",
        # "syn.hash.self_count",
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
    del anchor_hist["hash.self_count"]
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
        "syn.in_pct": synteny_pct,
        "syn.hashes.disambig": n_disambiguated,
        "syn.hashes.nonambig": n_nonambig,
        "syn.hashes.unassigned": n_assigned - in_synteny,
        "syn.hashes.unassigned_pct": unassigned_pct,
        "syn.fom": avg_ortho * 100.0 / n_proteomes,
    }
    return synteny_stats


def join_synteny_to_clusters(args, cluster_parent=None, mailbox_reader=None):
    """Read homology info from mailbox and join it to proteome file."""
    idx = args[0]
    cluster_path = cluster_parent / f"{idx}.parq"
    cluster = pd.read_parquet(cluster_path)
    n_cluster = len(cluster)
    with mailbox_reader(idx) as fh:
        synteny_frame = pd.read_csv(fh, sep="\t", index_col=0).convert_dtypes()
        in_synteny = len(synteny_frame)
    # cluster files are unusual in that I don't bother to version them,
    # so overlapping info has to be deleted each time
    clust_syn = concat_without_overlap(cluster, synteny_frame)
    write_tsv_or_parquet(clust_syn, cluster_path)
    anchor_count = clust_syn["syn.anchor_id"].value_counts()
    anchor_frag_counts = [0]
    for unused_id_tuple, subframe in clust_syn.groupby(
        by=["syn.anchor_id", "path"]
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


def dagchainer_id_to_int(ident):
    """Accept DAGchainer ids such as "cl1" and returns an integer."""
    if not ident.startswith("cl"):
        raise ValueError(f"Invalid ID {ident}.")
    id_val = ident[2:]
    if not id_val.isnumeric():
        raise ValueError(f"Non-numeric ID value in {ident}.")
    return int(id_val)
