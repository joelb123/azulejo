# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import fcntl
import json
import os
import statistics
import sys

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import click
import dask.bag as db
import gffpandas.gffpandas as gffpd
import numpy as np
import pandas as pd
import sh
from dask.diagnostics import ProgressBar
from loguru import logger

# module imports
from . import cli
from . import click_loguru
from .common import dotpath_to_path
from .common import parse_cluster_fasta
from .core import cleanup_fasta
from .core import usearch_cluster
from .inputs import TaxonomicInputTable
from .inputs import filepath_from_url
from .inputs import read_from_url
from .mailboxes import DataMailboxes

# global constants
PROTEOMES_FILE = "proteomes.tsv"
PROTEINS_FILE = "proteins.tsv"
HOMOLOGY_FILE = "proteins+homology.tsv"
SYNTENY_ENDING = "synteny.tsv"
PROXY_ENDING = "proxy.tsv"


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
        ending = HOMOLOGY_FILE
        file_type = "homology"
    else:
        ending = f"-{synteny}-{SYNTENY_ENDING}"
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


@cli.command()
@click_loguru.init_logger()
@click.option(
    "-s/-l",
    "--shorten_source/--no-shorten_source",
    default=True,
    is_flag=True,
    show_default=True,
    help="Remove invariant dotpaths in source IDs.",
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    default=True,
    show_default=True,
    help="Process in parallel.",
)
@click.argument("input_toml")
def ingest_data(shorten_source, input_toml, parallel):
    """Marshal protein and genome sequence information.


    IDs must correspond between GFF and FASTA files and must be unique across
    the entire set.
    """
    options = click_loguru.get_global_options()
    input_obj = TaxonomicInputTable(Path(input_toml), write_table=False)
    input_table = input_obj.input_table
    set_path = Path(input_obj.setname)
    arg_list = []
    for i, row in input_table.iterrows():
        arg_list.append((row["path"], row["fasta_url"], row["gff_url"],))
    if parallel:
        bag = db.from_sequence(arg_list)
    else:
        file_stats = []
    if not options.quiet:
        logger.info(f"Extracting FASTA/GFF info for {len(arg_list)} genomes:")
        ProgressBar().register()
    if parallel:
        file_stats = bag.map(
            read_fasta_and_gff, shorten_source=shorten_source
        ).compute()
    else:
        for args in arg_list:
            file_stats.append(
                read_fasta_and_gff(args, shorten_source=shorten_source)
            )
    del arg_list
    seq_frame = pd.DataFrame.from_dict([s[0] for s in file_stats]).set_index(
        "path"
    )
    frag_frame = pd.DataFrame.from_dict([s[1] for s in file_stats]).set_index(
        "path"
    )
    proteomes = pd.concat(
        [input_table.set_index("path"), frag_frame, seq_frame], axis=1
    )
    proteomes.drop(["fasta_url", "gff_url"], axis=1, inplace=True)
    proteomes = sort_proteome_frame(proteomes)
    proteome_table_path = set_path / PROTEOMES_FILE
    logger.info(
        f'Writing table of proteomes to "{proteome_table_path}", edit it to change preferences'
    )
    logger.info(
        "This is also the time to put DNA fragments on a common naming basis."
    )
    proteomes.to_csv(proteome_table_path, sep="\t")


def sort_proteome_frame(df):
    """Sort a proteome frame by preference and frag.max and renumber."""
    if df.index.name == "path":
        df["path"] = df.index
    df.sort_values(
        by=["preference", "frag.max"], ascending=[True, False], inplace=True
    )
    df["order"] = range(len(df))
    df.set_index("order", inplace=True)
    return df


def read_fasta_and_gff(args, shorten_source=True):
    """Read corresponding sequence and position files and construct consolidated tables."""
    dotpath, fasta_url, gff_url = args
    out_path = dotpath_to_path(dotpath)

    with read_from_url(fasta_url) as fasta_fh:
        unused_stem, unused_path, prop_frame, file_stats = cleanup_fasta(
            out_path, fasta_fh, dotpath, write_fasta=False, write_stats=False
        )
    # logger.debug(f"Reading GFF file {gff_url}.")
    with filepath_from_url(gff_url) as local_gff_file:
        annotation = gffpd.read_gff3(local_gff_file)
    mrnas = annotation.filter_feature_of_type(["mRNA"]).attributes_to_columns()
    del annotation
    mrnas.drop(
        mrnas.columns.drop(
            ["seq_id", "start", "strand", "ID"]
        ),  # drop EXCEPT these
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
    # Drop any mrnas not found in sequence file, e.g., zero-length
    mrnas = mrnas[mrnas.index.isin(prop_frame.index)]
    # sort fragments by largest value
    frag_counts = mrnas["frag_id"].value_counts()
    frag_frame = pd.DataFrame()
    frag_frame["counts"] = frag_counts
    frag_frame["idx"] = range(len(frag_frame))
    frag_frame["frag_id"] = frag_frame.index
    frag_frame["new_name"] = ""
    frag_frame.set_index(["idx"], inplace=True)
    frag_stats = {
        "path": dotpath,
        "frag.n": len(frag_frame),
        "frag.max": frag_counts[0],
    }
    frag_count_path = out_path / "fragments.tsv"
    if not frag_count_path.exists():
        # logger.debug(f"Writing fragment stats to file {frag_count_path}")
        frag_frame.to_csv(frag_count_path, sep="\t")
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
    joined_path = out_path / PROTEINS_FILE
    mrnas = mrnas.join(prop_frame)
    # logger.debug(f"Writing protein info file {joined_path}.")
    mrnas.to_csv(joined_path, sep="\t")
    return file_stats, frag_stats


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
def annotate_homology(identity, setname, parallel):
    """Add homology cluster info."""
    options = click_loguru.get_global_options()
    set_path = Path(setname)
    file_stats_path = set_path / PROTEOMES_FILE
    proteomes = pd.read_csv(file_stats_path, index_col=0, sep="\t")
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
    n_clusters, run_stats, cluster_hist = usearch_cluster.callback(
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
        mb_dir_path=(set_path / "mailboxes"),
        delete_on_exit=False,
    )
    mb.write_headers("\tcluster\tadj_group\n")
    cluster_paths = [
        set_path / "homology" / f"{i}.fa" for i in range(n_clusters)
    ]
    if parallel:
        bag = db.from_sequence(cluster_paths)
    else:
        cluster_stats = []
    if not options.quiet:
        logger.info(
            f"Calculating MSAs and trees for {len(cluster_paths)} homology clusters:"
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
    cluster_frame = pd.DataFrame.from_dict(clusters_dict).transpose()
    del clusters_dict
    cluster_frame.sort_index(inplace=True)
    grouping_dict = {}
    for i in range(n_proteomes):  # keep numbering of single-file clusters
        grouping_dict[f"[{i}]"] = i
    grouping_dict[str(list(range(n_proteomes)))] = 0
    for n_members, subframe in cluster_frame.groupby(["n_memb"]):
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
        keyfile = f"groupkeys-{n_members}.tsv"
        logger.debug(
            f"Writing group keys for group size {n_members} to {keyfile}"
        )
        member_counts.to_csv(set_path / keyfile, sep="\t")
    cluster_frame["members"] = cluster_frame["members"].map(grouping_dict)
    cluster_frame = cluster_frame.rename(columns={"members": "group_key"})
    n_adj = cluster_frame["n_adj"].sum()
    adj_pct = n_adj * 100.0 / n_clust_genes
    n_adj_clust = sum(cluster_frame["adj_groups"] != 0)
    adj_clust_pct = n_adj_clust * 100.0 / len(cluster_frame)
    logger.info(
        f"{n_adj} ({adj_pct:.1f}%) out of {n_clust_genes}"
        + " clustered genes are adjacent"
    )
    logger.info(
        f"{n_adj_clust} ({adj_clust_pct:.1f}%) out of "
        + f"{len(cluster_frame)} clusters contain adjacency"
    )
    cluster_file_name = f"clusters.tsv"
    logger.debug(f"Writing cluster file to {cluster_file_name}")
    cluster_frame.to_csv(set_path / cluster_file_name, sep="\t")
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
    mb.__del__(force=True)
    homo_frame = pd.DataFrame.from_dict(homo_stats)
    homo_frame.set_index(["idx"], inplace=True)
    homo_frame.sort_index(inplace=True)
    logger.info("Homology cluster coverage:")
    with pd.option_context(
        "display.max_rows", None, "display.float_format", "{:,.2f}%".format
    ):
        logger.info(homo_frame)
    proteomes = pd.concat([proteomes, homo_frame], axis=1)
    proteomes.to_csv(set_path / PROTEOMES_FILE, sep="\t", float_format="%5.2f")


def write_concatenated_protein_fasta(args):
    row, concat_fasta_path = args
    """Read peptide sequences from info file and write them out."""
    dotpath = row["path"]
    phylogeny_dict = {"idx": row.name, "path": dotpath}
    for n in [name for name in row.index if name.startswith("phy.")]:
        phylogeny_dict[n] = row[n]
    inpath = dotpath_to_path(dotpath) / PROTEINS_FILE
    prot_info = pd.read_csv(inpath, index_col=0, sep="\t")
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
    cluster_frame = pd.DataFrame.from_dict(prop_dict).transpose()
    cluster_frame["idx"] = cluster_frame["path"].map(file_dict)
    cluster_frame.sort_values(by=["idx", "frag_id", "frag_pos"], inplace=True)
    cluster_frame["adj_group"] = ""
    adjacency_group = 0
    was_adj = False
    for unused_group_id, subframe in cluster_frame.groupby(
        by=["idx", "frag_id"]
    ):
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
                    cluster_frame.loc[last_ID, "adj_group"] = str(
                        adjacency_group
                    )
                was_adj = True
                cluster_frame.loc[gene_id, "adj_group"] = str(adjacency_group)
            else:
                if was_adj:
                    adjacency_group += 1
                    was_adj = False
            last_pos = row["frag_pos"]
            last_ID = gene_id
    if was_adj:
        adjacency_group += 1
    idx_values = cluster_frame["idx"].value_counts()
    idx_list = list(idx_values.index)
    idx_list.sort()
    cluster_frame.to_csv(outdir / f"{cluster_id}.tsv", sep="\t")
    n_adj = sum(cluster_frame["adj_group"] != "")
    cluster_dict = {
        "size": len(cluster_frame),
        "n_memb": len(idx_values),
        "members": str(idx_list),
        "n_adj": n_adj,
        "adj_groups": adjacency_group,
    }
    for group_id, subframe in cluster_frame.groupby(by=["idx"]):
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
    proteins = pd.read_csv(
        protein_parent / PROTEINS_FILE, index_col=0, sep="\t"
    )
    n_proteins = len(proteins)
    with mailbox_reader(idx) as fh:
        homology_frame = pd.read_csv(fh, index_col=0, sep="\t")
        clusters_in_proteome = len(homology_frame)
    # homology_frame["adj_group"] = int(homology_frame["adj_group"])
    # homology_frame["cluster"] = int(homology_frame["cluster"])
    homology_frame = homology_frame[["cluster", "adj_group"]]
    proteome_frame = pd.concat([proteins, homology_frame], axis=1)
    proteome_frame.to_csv(protein_parent / HOMOLOGY_FILE, sep="\t")
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
        infoobj = pd.read_csv(infofile, index_col=0, sep="\t")
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
                f"subcluster {subcluster} is not unique w.r.t. proteome {list(stems)[best_choice]}."
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

    prefs is an optional list of proteome stems in order of preference in the proxy calc.
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
    proxy_filename = f"{setname}-{synteny_type}-{PROXY_ENDING}"
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
        f"{setname}-{synteny_type}-downselected-{PROXY_ENDING}"
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
