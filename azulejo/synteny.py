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
from .common import MERGED_HASH_FILE
from .common import PROTEOMOLOGY_FILE
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
@click.argument("setname")
def synteny_anchors(k, peatmer, setname):
    """Calculate synteny anchors."""
    options = click_loguru.get_global_options()
    set_path = Path(setname)
    file_stats_path = set_path / PROTEOMOLOGY_FILE
    proteomes = pd.read_csv(file_stats_path, index_col=0, sep="\t")
    n_proteomes = len(proteomes)
    if peatmer:
        synteny_type = "p"
    else:
        synteny_type = "kr"
    logger.info(
        f"Calculating {synteny_type}-mer of length ${k} synteny blocks with "
        f" for {n_proteomes} proteomes"
    )
    hash_mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "hash_merge"),
    )
    hash_mb.write_headers("hash\n")
    for idx, row in proteomes.iterrows():
        dotpath = row["path"]
        outpath = dotpath_to_path(dotpath)
        h = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
        h["nan_group"] = ((h["cluster"].isnull()).astype(int).cumsum() + 1) * (
            ~h["cluster"].isnull()
        )
        hasher = SyntenyBlockHasher(k=k, peatmer=peatmer)
        hash_name = hasher.hash_name()
        synteny_frame_list = []
        logger.debug(f"Calculating synteny blocks using {hash_name} .")
        for unused_id_tuple, subframe in h.groupby(
            by=["frag_id", "nan_group"]
        ):
            synteny_frame_list.append(hasher.calculate(subframe["cluster"]))
        synteny_frame = pd.concat(synteny_frame_list, axis=0)
        del synteny_frame_list
        hash_series = synteny_frame[hash_name]
        hash_counts = hash_series.value_counts()
        synteny_frame["syn.selfcount"] = synteny_frame[hash_name].map(
            hash_counts
        )
        write_tsv_or_parquet(synteny_frame, outpath / SYNTENY_FILE)
        # Do histogram of hashes
        hash_hist = pd.DataFrame(hash_counts.value_counts()).sort_index()
        hash_hist["pct_hashes"] = (
            hash_hist[hash_name] * 100.0 / len(synteny_frame)
        )
        write_tsv_or_parquet(outpath / HASH_HIST_FILE)
        # Write out sorted list of hash values
        unique_hashes = hash_series.unique().to_numpy().astype("uint32")
        unique_hashes.sort()
        with hash_mb.locked_open_for_write(idx) as fh:
            np.savetxt(fh, unique_hashes, fmt="%d")
    logger.info("Merging hashes")
    merger = ExternalMerge(
        file_path_func=hash_mb.path_to_mailbox, n_merge=n_proteomes
    )
    merger.init("hash")
    merged_hashes = merger.merge()
    hash_mb.delete()
    write_tsv_or_parquet(merged_hashes, set_path / MERGED_HASH_FILE)
    logger.info("Putting merged hashes into proteomes")
    for idx, row in proteomes.iterrows():
        dotpath = row["path"]
        outpath = dotpath_to_path(dotpath)
        syn = read_tsv_or_parquet(outpath / SYNTENY_FILE)
        syn = syn.join(merged_hashes, on=hash_name)
        homology = read_tsv_or_parquet(outpath / HOMOLOGY_FILE)
        syn = pd.concat([homology, syn], axis=1)
        syn["syn.id"] = ""
        id_col_no = syn.columns.get_loc("syn.id")
        ortho_col_no = syn.columns.get_loc("syn.ortho")
        cluster_col_no = syn.columns.get_loc("cluster")
        for ortho, subframe in syn.groupby(by=["syn.ortho"]):
            for unused_i, row in subframe.iterrows():
                footprint = row["syn.footprint"]
                direction = row["syn.direction"]
                cluster = None
                extra_indexes = 0
                if direction == 1:
                    k_list = list(reversed(range(footprint)))
                else:
                    k_list = list(range(footprint))
                prefix = row["syn.prefix"]
                row_no = syn.index.get_loc(row.name)
                for offset in range(footprint):
                    if peatmer:
                        cur_cluster = syn.iloc[row_no + offset, cluster_col_no]
                        if cur_cluster != cluster:
                            cluster = cur_cluster
                            sub_id = k_list.pop()
                    else:
                        sub_id = k_list.pop()
                    syn.iloc[
                        row_no + offset, id_col_no
                    ] = f"{ortho}.{int(prefix)}"
                    syn.iloc[row_no + offset, ortho_col_no] = ortho
        del (
            syn[hash_name],
            syn["syn.footprint"],
            syn["cluster"],
            syn["syn.prefix"],
        )
        write_tsv_or_parquet(
            syn, outpath / SYNTENY_FILE,
        )
        in_synteny = (syn["syn.id"] != "").sum()
        ambig = (syn["syn.selfcount"] != 1).sum()
        synteny_pct = in_synteny * 100.0 / len(syn)
        unambig_pct = (in_synteny - ambig) * 100.0 / len(syn)
        synteny_stats = {
            "path": dotpath,
            "genes": len(syn),
            "syntenic": in_synteny,
            "synteny_ambig": ambig,
            "synteny_pct": synteny_pct,
            "unambig_pct": unambig_pct,
        }
        logger.info(synteny_stats)
    logger.info("done with synteny")


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
