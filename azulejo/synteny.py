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
from .common import HASH_HIST_FILE
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

    def hash_name(self):
        """Return the string name of the hash function."""
        if self.peatmer:
            return f"syn.hash.peatmer{self.k}"
        return f"syn.hash.kmer{self.k}"

    def calculate(self, cluster_series):
        """Return an array of synteny block hashes data."""
        # Maybe the best code I've ever written--JB
        ids = []
        directions = array.array("h")  # signed 16-byte
        hashes = array.array("L")  # unsigned 32-bit
        # Set up indirect indexing over cluster list
        vec = cluster_series.to_numpy().astype(int)
        if self.peatmer:
            unequal_idxes = np.append(
                np.where(vec[1:] != vec[:-1]), len(vec) - 1
            )
            runlengths = np.diff(np.append(-1, unequal_idxes))
            positions = np.cumsum(np.append(0, runlengths))[:-1]
            n_iter = len(positions) - self.k + 1
            footprints = pd.array(
                [runlengths[i : i + self.k].sum() for i in range(n_iter)],
                dtype=pd.UInt32Dtype(),
            )
        else:
            n_elements = len(cluster_series)
            n_iter = n_elements - self.k + 1
            positions = np.arange(n_elements)
            footprints = pd.array([self.k] * (n_iter), dtype=pd.UInt32Dtype())
        for start in range(n_iter):
            kmer = vec[positions[start : start + self.k]]
            fwd_hash = xxhash.xxh32_intdigest(kmer.tobytes())
            rev_hash = xxhash.xxh32_intdigest(np.flip(kmer.tobytes()))
            if fwd_hash <= rev_hash:
                hashes.append(fwd_hash)
                directions.append(1)
            else:
                hashes.append(rev_hash)
                directions.append(-1)
            ids.append(cluster_series.index[start])
        hash_arr = pd.array(hashes, dtype=pd.UInt32Dtype())
        directions_ser = pd.Categorical(directions)
        return pd.DataFrame(
            [directions_ser, footprints, hash_arr],
            columns=["syn.direction", "syn.footprint", self.hash_name(),],
            index=ids,
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
    hash_name = hasher.hash_name()
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
            f"Calculating synteny anchors using the {hash_name} function"
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
    merged_hashes = merger.merge()
    hash_mb.delete()
    ret_list = []
    if not options.quiet:
        logger.info(
            f"Merging {len(merged_hashes)} synteny anchors into {n_proteomes} proteomes"
        )
        ProgressBar().register()
    if parallel:
        ret_list = bag.map(
            merge_synteny_hashes,
            merged_hashes=merged_hashes,
            hash_name=hash_name,
        ).compute()
    else:
        for args in arg_list:
            ret_list.append(
                merge_synteny_hashes(
                    args, merged_hashes=merged_hashes, hash_name=hash_name
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
    h = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
    h["nan_group"] = ((h["cluster"].isnull()).astype(int).cumsum() + 1) * (
        ~h["cluster"].isnull()
    )
    hash_name = hasher.hash_name()
    synteny_frame_list = []
    for unused_id_tuple, subframe in h.groupby(by=["frag.id", "nan_group"]):
        synteny_frame_list.append(hasher.calculate(subframe["cluster"]))
    synteny_frame = pd.concat(synteny_frame_list, axis=0)
    del synteny_frame_list
    hash_series = synteny_frame[hash_name]
    hash_counts = hash_series.value_counts()
    synteny_frame["syn.selfcount"] = synteny_frame[hash_name].map(hash_counts)
    write_tsv_or_parquet(synteny_frame, outpath / SYNTENY_FILE)
    # Do histogram of hashes
    hash_hist = pd.DataFrame(hash_counts.value_counts()).sort_index()
    hash_hist["pct_hashes"] = hash_hist[hash_name] * 100.0 / len(synteny_frame)
    write_tsv_or_parquet(hash_hist, outpath / HASH_HIST_FILE)
    # Write out sorted list of hash values
    unique_hashes = hash_series.unique().to_numpy().astype("uint32")
    unique_hashes.sort()
    with mailboxes.locked_open_for_write(idx) as fh:
        np.savetxt(fh, unique_hashes, fmt="%d")
    return len(unique_hashes)


def merge_synteny_hashes(args, merged_hashes=None, hash_name=None):
    """Merge synteny hashes into proteomes."""
    idx, dotpath = args
    outpath = dotpath_to_path(dotpath)
    syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
    syn = syn.join(merged_hashes, on=hash_name)
    homology = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
    syn = pd.concat([homology, syn], axis=1)
    syn["i"] = range(len(syn))
    synteny_id_arr = pd.array([pd.NA] * len(syn), dtype=pd.Int32Dtype())
    ortho_arr = pd.array([pd.NA] * len(syn), dtype=pd.Int32Dtype())
    for ortho, subframe in syn.groupby(by=["syn.ortho"]):
        for unused_i, row in subframe.iterrows():
            footprint = row["syn.footprint"]
            prefix = row["syn.prefix"]
            row_no = row["i"]
            synteny_id_arr[row_no : row_no + footprint] = prefix
            ortho_arr[row_no : row_no + footprint] = ortho
    syn["syn.id"] = synteny_id_arr
    syn["syn.ortho"] = ortho_arr
    del (
        syn[hash_name],
        syn["i"],
        syn["syn.prefix"],
    )
    write_tsv_or_parquet(
        syn, outpath / SYNTENY_FILE,
    )
    n_genes = len(syn)
    in_synteny = n_genes - syn["syn.id"].isnull().sum()
    n_assigned = n_genes - syn["cluster"].isnull().sum()
    ambig = (syn["syn.selfcount"] != 1).sum()
    synteny_pct = in_synteny * 100.0 / n_assigned
    unambig_pct = (in_synteny - ambig) * 100.0 / n_assigned
    synteny_stats = {
        "idx": idx,
        "path": dotpath,
        "hom.assign": n_assigned,
        "syn.assign": in_synteny,
        "synt.ambig": ambig,
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
        synteny_name = f"{stem}-{synteny_hash_name}{SYNTENY_ENDING}"
        logger.debug(
            f"Writing {synteny_hash_name} synteny frame {synteny_name}."
        )
        homology_frame.to_csv(set_path / synteny_name, sep="\t")
