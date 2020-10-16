# -*- coding: utf-8 -*-
"""Homology (sequence similarity) operations."""
# standard library imports
import fcntl
import json
import os
import shutil
import sys
from pathlib import Path

# third-party imports
import dask.bag as db
import pandas as pd
from dask.diagnostics import ProgressBar

# first-party imports
import sh

# module imports
from .common import CLUSTER_FILETYPE
from .common import CLUSTERS_FILE
from .common import EXTERNAL_CLUSTERS_FILE
from .common import FRAGMENTS_FILE
from .common import HOMOLOGY_FILE
from .common import PROTEINS_FILE
from .common import PROTEOMES_FILE
from .common import PROTEOMOLOGY_FILE
from .common import SPINNER_UPDATE_PERIOD
from .common import TrimmableMemoryMap
from .common import calculate_adjacency_group
from .common import dotpath_to_path
from .common import group_key_filename
from .common import logger
from .common import read_tsv_or_parquet
from .common import sort_proteome_frame
from .common import write_tsv_or_parquet
from .core import homology_cluster
from .mailboxes import DataMailboxes

# global constants
HOMOLOGY_COLS = ["hom.cluster", "hom.cl_size"]


def cluster_build_trees(
    identity, set_name, cluster_file=None, click_loguru=None
):
    """Calculate homology clusters, MSAs, trees."""
    options = click_loguru.get_global_options()
    user_options = click_loguru.get_user_global_options()
    parallel = user_options["parallel"]
    set_path = Path(set_name)
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
    frags["frag.idx"] = pd.array(frags.index, dtype=pd.UInt32Dtype())
    frag_frames = {}
    for dotpath, subframe in frags.groupby(by=["path"]):
        frag_frames[dotpath] = subframe.copy().set_index("frag.orig_id")
    arg_list = []
    concat_fasta_path = set_path / "proteins.fa"
    for i, row in proteomes.iterrows():
        arg_list.append((row, concat_fasta_path, frag_frames[row["path"]]))
    file_idx = {}
    stem_dict = {}
    for i, row in proteomes.iterrows():
        stem = row["path"]
        file_idx[stem] = i
        stem_dict[i] = stem
    if cluster_file is None:
        if concat_fasta_path.exists():
            concat_fasta_path.unlink()
        if not options.quiet:
            logger.info(
                f"Renaming fragments and concatenating sequences for {len(arg_list)}"
                " proteomes:"
            )
        for args in arg_list:
            write_protein_fasta(args)
        del arg_list
        cwd = Path.cwd()
        os.chdir(set_path)
        n_clusters, run_stats, cluster_hist = homology_cluster(
            "proteins.fa",
            identity,
            write_ids=True,
            delete=False,
            cluster_stats=False,
            outname="homology",
            click_loguru=click_loguru,
        )
        log_path = Path("homology.log")
        log_dir_path = Path("logs")
        log_dir_path.mkdir(exist_ok=True)
        shutil.copy2(log_path, "logs/homology.log")
        log_path.unlink()
        os.chdir(cwd)
        logger.info(f"Number of clusters: {n_clusters}")
        del cluster_hist
        del run_stats
        concat_fasta_path.unlink()
    else:  # use pre-existing clusters
        inclusts = pd.read_csv(cluster_file, sep="\t")
        for col in ["cluster_id", "members"]:
            if col not in inclusts.columns:
                logger.error(
                    f'Column named "{col}" not found in external homology cluster file'
                )
                sys.exit(1)
        cluster_counts = inclusts["cluster_id"].value_counts()
        cluster_map = pd.Series(
            range(len(cluster_counts)), index=cluster_counts.index
        )
        cluster_ids = inclusts["cluster_id"].map(cluster_map)
        cluster_sizes = inclusts["cluster_id"].map(cluster_counts)
        predef_clusters = pd.DataFrame(
            {
                "cluster_id": cluster_ids,
                "size": cluster_sizes,
                "members": inclusts["members"],
            }
        )
        predef_clusters.sort_values(by=["cluster_id"], inplace=True)
        predef_clusters.drop(
            predef_clusters[predef_clusters["size"] < 2].index,
            axis=0,
            inplace=True,
        )
        n_clusters = predef_clusters["cluster_id"].max() + 1
        predef_clusters.index = range(len(predef_clusters))
        external_cluster_path = set_path / EXTERNAL_CLUSTERS_FILE
        logger.info(
            f"Writing {external_cluster_path} with {len(predef_clusters)} genes"
            + f" in {n_clusters} homology clusters"
        )
        predef_clusters.to_csv(external_cluster_path, sep="\t")
        del cluster_counts, cluster_map, cluster_sizes, inclusts
        homology_path = set_path / "homology"
        homology_path.mkdir(exist_ok=True)
        if not options.quiet:
            logger.info(
                f"Creating cluster files for for {len(arg_list)}" " proteomes:"
            )
        proteome_no = 0
        for args in arg_list:
            logger.info(f"doing proteome {proteome_no}")
            write_protein_fasta(
                args, fasta_dir=homology_path, clusters=predef_clusters
            )
            proteome_no += 1
        del arg_list
        logger.info(
            "Checking that all cluster files are present (gene-id mismatch)"
        )
        missing_files = False
        for i in range(n_clusters):
            if not (homology_path / f"{i}.fa").exists():
                logger.error(f"External cluster {i} is missing.")
                missing_files = True
        if missing_files:
            sys.exit(1)
    #
    # Write homology info back into proteomes
    #
    click_loguru.elapsed_time("Alignment/tree-building")
    hom_mb = DataMailboxes(
        n_boxes=n_proteomes,
        mb_dir_path=(set_path / "mailboxes" / "clusters2proteomes"),
        file_extension="tsv",
    )
    hom_mb.write_tsv_headers(HOMOLOGY_COLS)
    cluster_paths = [
        set_path / "homology" / f"{i}.fa" for i in range(n_clusters)
    ]
    bag = db.from_sequence(cluster_paths)
    cluster_stats = []
    if not options.quiet:
        logger.info(
            f"Calculating MSAs and trees for {len(cluster_paths)} homology"
            " clusters:"
        )
        ProgressBar(dt=SPINNER_UPDATE_PERIOD).register()
    if parallel:
        cluster_stats = bag.map(
            parse_cluster,
            file_dict=file_idx,
            file_writer=hom_mb.locked_open_for_write,
        )
    else:
        for clust_fasta in cluster_paths:
            cluster_stats.append(
                parse_cluster(
                    clust_fasta,
                    file_dict=file_idx,
                    file_writer=hom_mb.locked_open_for_write,
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
        member_counts = pd.DataFrame(subframe["n_members"].value_counts())
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
    clusters["n_members"] = clusters["n_members"].map(grouping_dict)
    clusters = clusters.rename(columns={"n_members": "group_key"})
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
    click_loguru.elapsed_time("Joining")
    arg_list = []
    for i, row in proteomes.iterrows():
        arg_list.append((i, dotpath_to_path(row["path"]),))
    bag = db.from_sequence(arg_list)
    hom_stats = []
    if not options.quiet:
        logger.info(f"Joining homology info to {n_proteomes} proteomes:")
        ProgressBar(dt=SPINNER_UPDATE_PERIOD).register()
    if parallel:
        hom_stats = bag.map(
            join_homology_to_proteome, mailbox_reader=hom_mb.open_then_delete
        ).compute()
    else:
        for args in arg_list:
            hom_stats.append(
                join_homology_to_proteome(
                    args, mailbox_reader=hom_mb.open_then_delete
                )
            )
    hom_mb.delete()
    hom_frame = pd.DataFrame.from_dict(hom_stats)
    hom_frame.set_index(["prot.idx"], inplace=True)
    hom_frame.sort_index(inplace=True)
    logger.info("Homology cluster coverage:")
    with pd.option_context(
        "display.max_rows", None, "display.float_format", "{:,.2f}%".format
    ):
        logger.info(hom_frame)
    proteomes = pd.concat([proteomes, hom_frame], axis=1)
    write_tsv_or_parquet(
        proteomes, set_path / PROTEOMOLOGY_FILE, float_format="%5.2f"
    )
    click_loguru.elapsed_time(None)


def write_protein_fasta(args, clusters=None, fasta_dir=None):
    """Read peptide sequences from info file and write them out."""
    row, concat_fasta_path, frags = args
    dotpath = row["path"]
    phylogeny_dict = {"prot.idx": row.name, "path": dotpath}
    for phy_prop in [name for name in row.index if name.startswith("phy.")]:
        phylogeny_dict[phy_prop] = row[phy_prop]
    inpath = dotpath_to_path(dotpath)
    prot_info = read_tsv_or_parquet(inpath / PROTEINS_FILE)
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
    write_tsv_or_parquet(prot_info, inpath / HOMOLOGY_FILE)
    # include phylogeny info in per-sequence info
    for prop in phylogeny_dict:
        prot_info[prop] = phylogeny_dict[prop]
    # write concatenated sequence info
    if clusters is None:
        fasta_path = concat_fasta_path
        info_to_fasta(None, fasta_path, append=True, infoobj=prot_info)
    else:
        for cluster_id, subframe in clusters.groupby(by=["cluster_id"]):
            cluster_info = prot_info[prot_info.index.isin(subframe["members"])]
            fasta_path = fasta_dir / f"{cluster_id}.fa"
            info_to_fasta(None, fasta_path, append=True, infoobj=cluster_info)


def parse_cluster(
    fasta_path, file_dict=None, file_writer=None, neighbor_joining=False
):
    """Parse cluster FASTA headers to create cluster table.."""
    cluster_id = fasta_path.name[:-3]
    outdir = fasta_path.parent
    clusters = parse_cluster_fasta(fasta_path)
    if len(clusters) < 2:
        # fasta_path.unlink()
        logger.error(f"Singleton Cluster {cluster_id} is size {len(clusters)}")
        cluster_dict = {
            "size": len(clusters),
            "n_memb": None,
            "n_members": None,
            "n_adj": None,
            "adj_groups": None,
        }
        return int(cluster_id)
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
    # fasta_path.unlink()
    clusters["prot.idx"] = clusters["path"].map(file_dict)
    clusters.sort_values(by=["prot.idx", "frag.id", "frag.pos"], inplace=True)
    n_adj, adj_gr_count, unused_adj_group = calculate_adjacency_group(
        clusters["frag.pos"], clusters["frag.idx"]
    )
    idx_values = clusters["prot.idx"].value_counts()
    idx_list = list(idx_values.index)
    idx_list.sort()
    write_tsv_or_parquet(clusters, outdir / f"{cluster_id}.{CLUSTER_FILETYPE}")
    cluster_dict = {
        "size": len(clusters),
        "n_memb": len(idx_values),
        "n_members": str(idx_list),
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
        with file_writer(group_id) as file_handle:
            proteome_frame.to_csv(file_handle, header=False, sep="\t")
    return int(cluster_id), cluster_dict


def parse_cluster_fasta(filepath, trim_dict=True):
    """Return FASTA headers as a dictionary of properties."""
    next_pos = 0
    properties_dict = {}
    memory_map = TrimmableMemoryMap(filepath)
    with memory_map.map() as mem_map:
        size = memory_map.size
        next_pos = mem_map.find(b">", next_pos)
        while next_pos != -1 and next_pos < size:
            eol_pos = mem_map.find(b"\n", next_pos)
            if eol_pos == -1:
                break
            space_pos = mem_map.find(b" ", next_pos + 1, eol_pos)
            if space_pos == -1:
                raise ValueError(
                    f"Header format is bad in {filepath} header"
                    f" {len(properties_dict)+1}"
                )
            cluster_id = mem_map[next_pos + 1 : space_pos].decode("utf-8")
            payload = json.loads(mem_map[space_pos + 1 : eol_pos])
            properties_dict[cluster_id] = payload
            if trim_dict:
                size = memory_map.trim(space_pos, eol_pos)
            next_pos = mem_map.find(b">", space_pos)
    cluster = (
        pd.DataFrame.from_dict(properties_dict).transpose().convert_dtypes()
    )
    return cluster


def join_homology_to_proteome(args, mailbox_reader=None):
    """Read homology info from mailbox and join it to proteome file."""
    idx, protein_parent = args
    proteins = pd.read_parquet(protein_parent / HOMOLOGY_FILE)
    n_proteins = len(proteins)
    with mailbox_reader(idx) as file_handle:
        homology_frame = pd.read_csv(
            file_handle, sep="\t", index_col=0
        ).convert_dtypes()
        clusters_in_proteome = len(homology_frame)
    proteome_frame = pd.concat([proteins, homology_frame], axis=1)
    write_tsv_or_parquet(proteome_frame, protein_parent / HOMOLOGY_FILE)
    return {
        "prot.idx": idx,
        "hom.clusters": clusters_in_proteome,
        "hom.cluster_pct": clusters_in_proteome * 100.0 / n_proteins,
    }


def info_to_fasta(infofile, fastafile, append, infoobj=None):
    """Convert infofile to FASTA file."""
    if infoobj is None:
        infoobj = read_tsv_or_parquet(infofile)
    if append:
        filemode = "a+"
    else:
        filemode = "w"
    with Path(fastafile).open(filemode) as file_handle:
        fcntl.flock(file_handle, fcntl.LOCK_EX)
        logger.debug(f"Writing to {fastafile} with mode {filemode}.")
        seqs = infoobj["prot.seq"].copy()
        del infoobj["prot.seq"]
        for gene_id, row in infoobj.iterrows():
            file_handle.write(f">{gene_id} {row.to_json()}\n")
            file_handle.write(f"{seqs[gene_id]}\n")
        fcntl.flock(file_handle, fcntl.LOCK_UN)
