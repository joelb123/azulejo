# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import sys

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import attr
import click
import dask.bag as db
import numpy as np
import pandas as pd
from dask.diagnostics import ProgressBar

# first-party imports
import sh
import xxhash
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import ANCHOR_HIST_FILE
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
                f"Workaround for minor error in shingle at hash {hash}, base"
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
    hasher = SyntenyBlockHasher(k=k, peatmer=peatmer)
    hash_mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "hash_merge"),
    )
    hash_mb.write_headers("hash\n")
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
    merged_hashes = merger.merge(
        count_key="tmp.ortho_count", ordinal_key="tmp.base"
    )
    merged_hashes["tmp.ortho_count"] *= k
    hash_mb.delete()
    ret_list = []
    if not options.quiet:
        logger.info(
            f"Merging {len(merged_hashes)} synteny anchors into {n_proteomes}"
            " proteomes"
        )
        ProgressBar().register()
    if parallel:
        ret_list = bag.map(
            merge_synteny_hashes, merged_hashes=merged_hashes, hasher=hasher,
        ).compute()
    else:
        for args in arg_list:
            ret_list.append(
                merge_synteny_hashes(
                    args, merged_hashes=merged_hashes, hasher=hasher
                )
            )
    synteny_stats = pd.DataFrame.from_dict(ret_list)
    synteny_stats = synteny_stats.set_index("idx").sort_index()
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.float_format",
        "{:,.2f}%".format,
    ):
        logger.info(synteny_stats)
    del synteny_stats["path"]
    proteomes = pd.concat([proteomes, synteny_stats], axis=1)
    write_tsv_or_parquet(proteomes, set_path / PROTEOSYN_FILE)


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
    # Write out sorted list of hash values
    unique_hashes = syn[hash_name].unique().to_numpy()
    unique_hashes.sort()
    with mailboxes.locked_open_for_write(idx) as fh:
        np.savetxt(fh, unique_hashes, fmt="%d")
    return len(unique_hashes)


def merge_synteny_hashes(args, merged_hashes=None, hasher=None):
    """Merge synteny hashes into proteomes."""
    hash_name = hasher.hash_name()
    plain_hash_name = hasher.hash_name(no_prefix=True)
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = syn.join(merged_hashes, on=hash_name)
    homology = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
    syn = pd.concat([homology, syn], axis=1)
    print(f"doing {dotpath}")
    n_prots = len(syn)
    syn["tmp.i"] = range(len(syn))
    shingled_vars = {
        plain_hash_name: np.array([np.nan] * n_prots),
        "direction": np.array([""] * n_prots),
        "self_count": np.array([np.nan] * n_prots),
        "footprint": np.array([np.nan] * n_prots),
        "frag_count": np.array([np.nan] * n_prots),
        "ortho_count": np.array([np.nan] * n_prots),
    }
    anchor_blocks = np.array([np.nan] * n_prots)
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
    in_synteny = n_prots - syn["syn.anchor_id"].isnull().sum()
    # n_assigned = n_genes - syn["hom.cluster"].isnull().sum()
    n_assigned = 1
    # Do histogram of blocks
    # anchor_hist = pd.DataFrame(anchor_counts.value_counts()).sort_index()
    # anchor_hist = anchor_hist.rename(columns={"syn.anchor_id": "self_count"})
    # anchor_hist["pct_anchors"] = anchor_hist["self_count"] * anchor_hist.index * 100.0 / n_assigned
    # write_tsv_or_parquet(anchor_hist, outpath / ANCHOR_HIST_FILE)
    ambig = (syn["tmp.self_count"] != 1).sum()
    synteny_pct = in_synteny * 100.0 / n_assigned
    unambig_pct = (in_synteny - ambig) * 100.0 / n_assigned
    synteny_stats = {
        "idx": idx,
        "path": dotpath,
        "hom.assign": n_assigned,
        "syn.anchors": in_synteny,
        "syn.ambig": ambig,
        "syn.assgn_pct": synteny_pct,
        "syn.unamb_pct": unambig_pct,
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
