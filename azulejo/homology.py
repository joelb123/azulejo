# -*- coding: utf-8 -*-
"""Homology (sequence similarity) operations."""
# standard library imports
import fcntl
import json
import os
import shutil
from pathlib import Path

# third-party imports
import click
import dask.bag as db
import numpy as np
import pandas as pd
from dask.diagnostics import ProgressBar

# first-party imports
import sh
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import CLUSTER_FILETYPE
from .common import CLUSTERS_FILE
from .common import DIRECTIONAL_CATEGORY
from .common import FRAGMENTS_FILE
from .common import HOMOLOGY_FILE
from .common import PROTEINS_FILE
from .common import PROTEOMES_FILE
from .common import PROTEOMOLOGY_FILE
from .common import UNRENAMED_PROTEINS_FILE
from .common import TrimmableMemoryMap
from .common import dotpath_to_path
from .common import group_key_filename
from .common import read_tsv_or_parquet
from .common import sort_proteome_frame
from .common import write_tsv_or_parquet
from .core import homology_cluster
from .mailboxes import DataMailboxes

# global constants

HOMOLOGY_COLS = ["hom.cluster", "hom.cl_size"]


@cli.command()
@click_loguru.init_logger()
@click.option(
    "--identity",
    "-i",
    default=0.0,
    help="Minimum sequence ID (0-1). [default: lowest]",
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    default=True,
    show_default=True,
    help="Process in parallel.",
)
@click.argument("setname")
def do_homology(identity, setname, parallel):
    """Calculate homology clusters, MSAs, trees."""
    options = click_loguru.get_global_options()
    set_path = Path(setname)
    # read and possibly update proteomes
    proteomes_path = set_path / PROTEOMES_FILE
    proteomes_in = read_tsv_or_parquet(proteomes_path)
    proteomes = sort_proteome_frame(proteomes_in)
    if not proteomes_in.equals(proteomes):
        logger.info("proteomes sort order changed, writing new proteomes file")
        write_tsv_or_parquet(proteomes, proteomes_path)
    n_proteomes = len(proteomes)
    # read and update fragment ID's
    frags = read_tsv_or_parquet(set_path / FRAGMENTS_FILE)
    frags["frag.idx"] = frags.index
    frag_frames = {}
    for dotpath, subframe in frags.groupby(by=["path"]):
        frag_frames[dotpath] = subframe.copy().set_index("frag.orig_id")
    # write concatenated stats
    concat_fasta_path = set_path / "proteins.fa"
    if concat_fasta_path.exists():
        concat_fasta_path.unlink()
    arg_list = []
    for i, row in proteomes.iterrows():
        arg_list.append((row, concat_fasta_path, frag_frames[row["path"]]))
    if not options.quiet:
        logger.info(
            f"Renaming frags and concatenating sequences for {len(arg_list)}"
            " proteomes:"
        )
    for args in arg_list:
        write_concatenated_protein_fasta(args)
    del arg_list
    file_idx = {}
    stem_dict = {}
    for i, row in proteomes.iterrows():
        stem = row["path"]
        file_idx[stem] = i
        stem_dict[i] = stem
    logger.debug("Doing cluster calculation.")
    cwd = Path.cwd()
    os.chdir(set_path)
    n_clusters, run_stats, cluster_hist = homology_cluster.callback(
        "proteins.fa",
        identity,
        write_ids=True,
        delete=False,
        cluster_stats=False,
        outname="homology",
    )
    log_path = Path("homology.log")
    shutil.copy2(log_path, "logs/homology.log")
    log_path.unlink()
    os.chdir(cwd)
    # logger.info(f"Stats of {n_clusters} clusters:")
    # logger.info(run_stats)
    # logger.info(f"\nCluster size histogram ({n_proteomes} proteomes):")
    # with pd.option_context(
    #     "display.max_rows", None, "display.float_format", "{:,.2f}%".format
    # ):
    #    logger.info(cluster_hist)
    del cluster_hist
    del run_stats
    concat_fasta_path.unlink()
    mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "clusters2proteomes"),
        file_extension="tsv",
    )
    mb.write_tsv_headers(HOMOLOGY_COLS)
    cluster_paths = [
        set_path / "homology" / f"{i}.fa" for i in range(n_clusters)
    ]
    if parallel:
        bag = db.from_sequence(cluster_paths)
    else:
        cluster_stats = []
    if not options.quiet:
        logger.info(
            f"Calculating MSAs and trees for {len(cluster_paths)} homology"
            " clusters:"
        )
        ProgressBar().register()
    if parallel:
        cluster_stats = bag.map(
            parse_cluster,
            file_dict=file_idx,
            file_writer=mb.locked_open_for_write,
        )
    else:
        for clust_fasta in cluster_paths:
            cluster_stats.append(
                parse_cluster(
                    clust_fasta,
                    file_dict=file_idx,
                    file_writer=mb.locked_open_for_write,
                )
            )
    n_clust_genes = 0
    clusters_dict = {}
    for cluster_id, cluster_dict in cluster_stats:
        n_clust_genes += cluster_dict["size"]
        clusters_dict[cluster_id] = cluster_dict
    del cluster_stats
    clusters = pd.DataFrame.from_dict(clusters_dict).transpose()
    del clusters_dict
    clusters.sort_index(inplace=True)
    grouping_dict = {}
    for i in range(n_proteomes):  # keep numbering of single-file clusters
        grouping_dict[f"[{i}]"] = i
    grouping_dict[str(list(range(n_proteomes)))] = 0
    for n_members, subframe in clusters.groupby(["n_memb"]):
        if n_members == 1:
            continue
        if n_members == n_proteomes:
            continue
        member_counts = pd.DataFrame(subframe["members"].value_counts())
        member_counts["key"] = range(len(member_counts))
        for newcol in range(n_members):
            member_counts[f"memb{newcol}"] = ""
        for member_string, row in member_counts.iterrows():
            grouping_dict[member_string] = row["key"]
            member_list = json.loads(member_string)
            for col in range(n_members):
                member_counts.loc[member_string, f"memb{col}"] = stem_dict[
                    member_list[col]
                ]
        member_counts = member_counts.set_index("key")
        write_tsv_or_parquet(
            member_counts, set_path / group_key_filename(n_members)
        )
    clusters["members"] = clusters["members"].map(grouping_dict)
    clusters = clusters.rename(columns={"members": "group_key"})
    n_adj = clusters["n_adj"].sum()
    adj_pct = n_adj * 100.0 / n_clust_genes
    n_adj_clust = sum(clusters["adj_groups"] != 0)
    adj_clust_pct = n_adj_clust * 100.0 / len(clusters)
    logger.info(
        f"{n_adj} ({adj_pct:.1f}%) out of {n_clust_genes}"
        + " clustered genes are adjacent"
    )
    logger.info(
        f"{n_adj_clust} ({adj_clust_pct:.1f}%) out of "
        + f"{len(clusters)} clusters contain adjacency"
    )
    write_tsv_or_parquet(clusters, set_path / CLUSTERS_FILE)
    # join homology cluster info to proteome info
    arg_list = []
    for i, row in proteomes.iterrows():
        arg_list.append((i, dotpath_to_path(row["path"]),))
    if parallel:
        bag = db.from_sequence(arg_list)
    else:
        homo_stats = []
    if not options.quiet:
        logger.info(f"Joining homology info to {n_proteomes} proteomes:")
        ProgressBar().register()
    if parallel:
        homo_stats = bag.map(
            join_homology_to_proteome, mailbox_reader=mb.open_then_delete
        ).compute()
    else:
        for args in arg_list:
            homo_stats.append(
                join_homology_to_proteome(
                    args, mailbox_reader=mb.open_then_delete
                )
            )
    mb.delete()
    homo_frame = pd.DataFrame.from_dict(homo_stats)
    homo_frame.set_index(["prot.idx"], inplace=True)
    homo_frame.sort_index(inplace=True)
    logger.info("Homology cluster coverage:")
    with pd.option_context(
        "display.max_rows", None, "display.float_format", "{:,.2f}%".format
    ):
        logger.info(homo_frame)
    proteomes = pd.concat([proteomes, homo_frame], axis=1)
    write_tsv_or_parquet(
        proteomes, set_path / PROTEOMOLOGY_FILE, float_format="%5.2f"
    )


def write_concatenated_protein_fasta(args):
    row, concat_fasta_path, frags = args
    """Read peptide sequences from info file and write them out."""
    dotpath = row["path"]
    phylogeny_dict = {"prot.idx": row.name, "path": dotpath}
    for n in [name for name in row.index if name.startswith("phy.")]:
        phylogeny_dict[n] = row[n]
    inpath = dotpath_to_path(dotpath)
    prot_info = read_tsv_or_parquet(inpath / UNRENAMED_PROTEINS_FILE)
    prot_info["frag.idx"] = prot_info["frag.id"].map(
        lambda oid: frags.loc[oid]["frag.idx"]
    )
    prot_info["frag.is_plas"] = prot_info["frag.id"].map(
        lambda oid: frags.loc[oid]["frag.is_plas"]
    )
    prot_info["frag.is_scaf"] = prot_info["frag.id"].map(
        lambda oid: frags.loc[oid]["frag.is_scaf"]
    )
    prot_info["frag.is_chr"] = prot_info["frag.id"].map(
        lambda oid: frags.loc[oid]["frag.is_chr"]
    )
    prot_info["frag.id"] = prot_info["frag.id"].map(
        lambda oid: frags.loc[oid]["frag.id"]
    )
    # Write out updated protein info
    write_tsv_or_parquet(prot_info, inpath / PROTEINS_FILE)
    # include phylogeny info in per-sequence info
    for prop in phylogeny_dict:
        prot_info[prop] = phylogeny_dict[prop]
    # write concatenated sequence info
    info_to_fasta.callback(
        None, concat_fasta_path, append=True, infoobj=prot_info
    )


def parse_cluster(
    fasta_path, file_dict=None, file_writer=None, neighbor_joining=False
):
    """Parse cluster FASTA headers to create cluster table.."""
    cluster_id = fasta_path.name[:-3]
    outdir = fasta_path.parent
    clusters = parse_cluster_fasta(fasta_path)
    if len(clusters) < 2:
        fasta_path.unlink()
        raise ValueError("Singleton Cluster")
    # calculate MSA and return guide tree
    muscle_args = [
        "-in",
        f"{outdir}/{cluster_id}.fa",
        "-out",
        f"{outdir}/{cluster_id}.faa",
        "-diags",
        "-sv",
        "-maxiters",
        "2",
        "-quiet",
        "-distance1",
        "kmer20_4",
    ]
    if len(clusters) >= 4:
        muscle_args += [
            "-tree2",
            f"{outdir}/{cluster_id}.nwk",
        ]
        if neighbor_joining:
            muscle_args += ["-cluster2", "neighborjoining"]  # adds 20%
    muscle = sh.Command("muscle")
    muscle(muscle_args)
    fasta_path.unlink()
    clusters["prot.idx"] = clusters["path"].map(file_dict)
    clusters.sort_values(by=["prot.idx", "frag.id", "frag.pos"], inplace=True)
    clusters["tmp.i"] = range(len(clusters))
    adj_gr_count = 0
    was_adj = False
    adj_group = np.array([np.nan] * len(clusters))
    for unused_group_id, subframe in clusters.groupby(
        by=["prot.idx", "frag.id"]
    ):
        if len(subframe) == 1:
            continue
        last_pos = -2
        last_row = None
        if was_adj:
            adj_gr_count += 1
        was_adj = False
        for unused_gene_id, row in subframe.iterrows():
            row_no = row["tmp.i"]
            if row["frag.pos"] == last_pos + 1:
                if not was_adj:
                    adj_group[last_row] = adj_gr_count
                was_adj = True
                adj_group[row_no] = adj_gr_count
            else:
                if was_adj:
                    adj_gr_count += 1
                    was_adj = False
            last_pos = row["frag.pos"]
            last_row = row_no
    if was_adj:
        adj_gr_count += 1
    clusters["adj_group"] = pd.array(adj_group, dtype=pd.UInt32Dtype())
    idx_values = clusters["prot.idx"].value_counts()
    idx_list = list(idx_values.index)
    idx_list.sort()
    write_tsv_or_parquet(clusters, outdir / f"{cluster_id}.{CLUSTER_FILETYPE}")
    n_adj = len(clusters) - clusters["adj_group"].isnull().sum()
    cluster_dict = {
        "size": len(clusters),
        "n_memb": len(idx_values),
        "members": str(idx_list),
        "n_adj": n_adj,
        "adj_groups": adj_gr_count,
    }
    for group_id, subframe in clusters.groupby(by=["prot.idx"]):
        proteome_frame = subframe.copy()
        proteome_frame["hom.cluster"] = cluster_id
        proteome_frame["hom.cl_size"] = len(idx_values)
        proteome_frame.drop(
            proteome_frame.columns.drop(HOMOLOGY_COLS),  # drop EXCEPT these
            axis=1,
            inplace=True,
        )
        with file_writer(group_id) as fh:
            proteome_frame.to_csv(fh, header=False, sep="\t")
    return int(cluster_id), cluster_dict


def parse_cluster_fasta(filepath, trim_dict=True):
    """Return FASTA headers as a dictionary of properties."""
    next_pos = 0
    properties_dict = {}
    memory_map = TrimmableMemoryMap(filepath)
    with memory_map.map() as mm:
        size = memory_map.size
        next_pos = mm.find(b">", next_pos)
        while next_pos != -1 and next_pos < size:
            eol_pos = mm.find(b"\n", next_pos)
            if eol_pos == -1:
                break
            space_pos = mm.find(b" ", next_pos + 1, eol_pos)
            if space_pos == -1:
                raise ValueError(
                    f"Header format is bad in {filepath} header"
                    f" {len(properties_dict)+1}"
                )
            id = mm[next_pos + 1 : space_pos].decode("utf-8")
            payload = json.loads(mm[space_pos + 1 : eol_pos])
            properties_dict[id] = payload
            if trim_dict:
                size = memory_map.trim(space_pos, eol_pos)
            next_pos = mm.find(b">", space_pos)
    cluster = (
        pd.DataFrame.from_dict(properties_dict).transpose().convert_dtypes()
    )
    cluster["frag.direction"] = pd.Categorical(
        cluster["frag.direction"], dtype=DIRECTIONAL_CATEGORY
    )
    for int_key in ["frag.start", "frag.pos", "prot.len", "prot.idx"]:
        cluster[int_key] = pd.array(cluster[int_key], dtype=pd.UInt32Dtype())
    return cluster


def join_homology_to_proteome(args, mailbox_reader=None):
    """Read homology info from mailbox and join it to proteome file."""
    idx, protein_parent = args
    proteins = pd.read_parquet(protein_parent / PROTEINS_FILE)
    n_proteins = len(proteins)
    with mailbox_reader(idx) as fh:
        homology_frame = pd.read_csv(
            fh, sep="\t", index_col=0
        ).convert_dtypes()
        clusters_in_proteome = len(homology_frame)
    proteome_frame = pd.concat([proteins, homology_frame], axis=1)
    write_tsv_or_parquet(proteome_frame, protein_parent / HOMOLOGY_FILE)
    return {
        "prot.idx": idx,
        "hom.clusters": clusters_in_proteome,
        "hom.cluster_pct": clusters_in_proteome * 100.0 / n_proteins,
    }


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--append/--no-append",
    "-a/-x",
    is_flag=True,
    default=True,
    help="Append to FASTA file.",
    show_default=True,
)
@click.argument("infofile")
@click.argument("fastafile")
def info_to_fasta(infofile, fastafile, append, infoobj=None):
    """Convert infofile to FASTA file."""
    if infoobj is None:
        infoobj = read_tsv_or_parquet(infofile)
    if append:
        filemode = "a+"
    else:
        filemode = "w"
    with Path(fastafile).open(filemode) as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        logger.debug(f"Writing to {fastafile} with mode {filemode}.")
        seqs = infoobj["prot.seq"].copy()
        del infoobj["prot.seq"]
        for gene_id, row in infoobj.iterrows():
            fh.write(f">{gene_id} {row.to_json()}\n")
            fh.write(f"{seqs[gene_id]}\n")
        fcntl.flock(fh, fcntl.LOCK_UN)
