# -*- coding: utf-8 -*-
"""Homology (sequence similarity) operations."""
# standard library imports
import fcntl
import json
import os
from pathlib import Path

# third-party imports
import click
import dask.bag as db
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
from .common import HOMOLOGY_FILE
from .common import PROTEINS_FILE
from .common import PROTEOMES_FILE
from .common import PROTEOMOLOGY_FILE
from .common import dotpath_to_path
from .common import group_key_filename
from .common import parse_cluster_fasta
from .common import read_tsv_or_parquet
from .common import write_tsv_or_parquet
from .core import homology_cluster
from .mailboxes import DataMailboxes


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
    file_stats_path = set_path / PROTEOMES_FILE
    proteomes = read_tsv_or_parquet(file_stats_path)
    n_proteomes = len(proteomes)
    # write concatenated stats
    concat_fasta_path = set_path / "proteins.fa"
    if concat_fasta_path.exists():
        concat_fasta_path.unlink()
    arg_list = []
    for i, row in proteomes.iterrows():
        arg_list.append((row, concat_fasta_path,))
    if not options.quiet:
        logger.info(f"Concatenating sequences for {len(arg_list)} proteomes:")
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
    os.chdir(cwd)
    logger.info(f"Stats of {n_clusters} clusters:")
    logger.info(run_stats)
    logger.info(f"\nCluster size histogram ({n_proteomes} proteomes):")
    with pd.option_context(
        "display.max_rows", None, "display.float_format", "{:,.2f}%".format
    ):
        logger.info(cluster_hist)
    del cluster_hist
    del run_stats
    concat_fasta_path.unlink()
    mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "clusters2proteomes"),
        file_extension="tsv",
    )
    mb.write_headers("\tadj_group\tcluster\n")
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
    homo_frame.set_index(["idx"], inplace=True)
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
    row, concat_fasta_path = args
    """Read peptide sequences from info file and write them out."""
    dotpath = row["path"]
    phylogeny_dict = {"idx": row.name, "path": dotpath}
    for n in [name for name in row.index if name.startswith("phy.")]:
        phylogeny_dict[n] = row[n]
    inpath = dotpath_to_path(dotpath) / PROTEINS_FILE
    prot_info = pd.read_parquet(inpath)
    for prop in phylogeny_dict:
        prot_info[prop] = phylogeny_dict[prop]
    info_to_fasta.callback(
        None, concat_fasta_path, append=True, infoobj=prot_info
    )


def parse_cluster(fasta_path, file_dict=None, file_writer=None):
    """Parse cluster FASTA headers to create cluster table.."""
    cluster_id = fasta_path.name[:-3]
    outdir = fasta_path.parent
    prop_dict = parse_cluster_fasta(fasta_path)
    if len(prop_dict) < 2:
        logger.warning(f"singleton cluster {fasta_path} removed")
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
    if len(prop_dict) >= 4:
        muscle_args += [
            "-tree2",
            f"{outdir}/{cluster_id}.nwk",
        ]  # ,  "-cluster2", "neighborjoining"] #adds 20%
    muscle = sh.Command("muscle")
    muscle(muscle_args)
    clusters = pd.DataFrame.from_dict(prop_dict).transpose()
    clusters["idx"] = clusters["path"].map(file_dict)
    clusters.sort_values(by=["idx", "frag_id", "frag_pos"], inplace=True)
    clusters["adj_group"] = ""
    adjacency_group = 0
    was_adj = False
    for unused_group_id, subframe in clusters.groupby(by=["idx", "frag_id"]):
        if len(subframe) == 1:
            continue
        last_pos = -2
        last_ID = None
        if was_adj:
            adjacency_group += 1
        was_adj = False
        for gene_id, row in subframe.iterrows():
            if row["frag_pos"] == last_pos + 1:
                if not was_adj:
                    clusters.loc[last_ID, "adj_group"] = str(adjacency_group)
                was_adj = True
                clusters.loc[gene_id, "adj_group"] = str(adjacency_group)
            else:
                if was_adj:
                    adjacency_group += 1
                    was_adj = False
            last_pos = row["frag_pos"]
            last_ID = gene_id
    if was_adj:
        adjacency_group += 1
    idx_values = clusters["idx"].value_counts()
    idx_list = list(idx_values.index)
    idx_list.sort()
    write_tsv_or_parquet(clusters, outdir / f"{cluster_id}.{CLUSTER_FILETYPE}")
    n_adj = sum(clusters["adj_group"] != "")
    cluster_dict = {
        "size": len(clusters),
        "n_memb": len(idx_values),
        "members": str(idx_list),
        "n_adj": n_adj,
        "adj_groups": adjacency_group,
    }
    for group_id, subframe in clusters.groupby(by=["idx"]):
        proteome_frame = subframe.copy()
        proteome_frame["cluster"] = cluster_id
        proteome_frame.drop(
            proteome_frame.columns.drop(
                ["adj_group", "cluster"]
            ),  # drop EXCEPT these
            axis=1,
            inplace=True,
        )
        with file_writer(group_id) as fh:
            proteome_frame.to_csv(fh, header=False, sep="\t")
    return int(cluster_id), cluster_dict


def join_homology_to_proteome(args, mailbox_reader=None):
    """Read homology info from mailbox and join it to proteome file."""
    idx, protein_parent = args
    proteins = pd.read_parquet(protein_parent / PROTEINS_FILE)
    n_proteins = len(proteins)
    with mailbox_reader(idx) as fh:
        homology_frame = read_tsv_or_parquet(fh)
        clusters_in_proteome = len(homology_frame)
    homology_frame = homology_frame[["cluster", "adj_group"]]
    proteome_frame = pd.concat([proteins, homology_frame], axis=1)
    write_tsv_or_parquet(proteome_frame, protein_parent / HOMOLOGY_FILE)
    return {
        "idx": idx,
        "clustered": clusters_in_proteome,
        "cluster_pct": clusters_in_proteome * 100.0 / n_proteins,
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
        for gene_id, row in infoobj.iterrows():
            row_dict = row.to_dict()
            seq = row_dict["prot.seq"]
            del row_dict["prot.seq"]
            json_row = json.dumps(row_dict, separators=(",", ":"))
            fh.write(f">{gene_id} {json_row}\n")
            fh.write(f"{seq}\n")
        fcntl.flock(fh, fcntl.LOCK_UN)
