# -*- coding: utf-8 -*-
#
# standard library imports
#
import os
import sys

#
# third-party imports
#
import click
import dask.bag as db
from dask.diagnostics import ProgressBar
import gffpandas.gffpandas as gffpd
import numpy as np
import pandas as pd
from Bio import SeqIO
from loguru import logger

#
# package imports
#
from . import cli
from . import click_loguru
from .common import *
from .core import usearch_cluster
from .core import cluster_set_name

# global constants
GFF_EXT = "gff3"
FASTA_EXT = "faa"
SEQ_ID_POS = 0
START_POS = 1
STRAND_POS = 2
PROTEIN_LEN_POS = 3
CLUSTER_POS = 4
FOOTPRINT_POS = 5
HASHDIR_POS = 6
SYNTENY_POS = 7
HOMOLOGY_ENDING = "-homology.tsv"
FILES_ENDING = "-files.tsv"
SYNTENY_ENDING = "-synteny.tsv"
PROXY_ENDING = "-proxy.tsv"

def kmer_synteny_block_func(k, frame):
    """Return a synteny block closure and its name."""
    frame_len = len(frame)
    frame_iter = iter(range(frame_len))
    null_return = (0, 0, 0, )
    def kmer_synteny():
        """Mappable function from index to (footprint, direction, hash) tuple."""
        first_index = next(frame_iter)
        cluster_list = []
        initial_seq_id = frame.iloc[first_index, SEQ_ID_POS]
        for idx in range(first_index, first_index+k):
            if idx+1 > frame_len:
                return null_return
            cluster =  frame.iloc[idx, CLUSTER_POS]
            if frame.iloc[idx, SEQ_ID_POS] != initial_seq_id or cluster == -1:
                return null_return
            cluster_list.append(cluster)
        fwd_hash = hash(tuple(cluster_list))
        rev_hash = hash(tuple(reversed(cluster_list)))
        if fwd_hash > rev_hash:
            hashval = fwd_hash
            direction = 1
        else:
            hashval = rev_hash
            direction = -1
        return k, direction, hashval

    return kmer_synteny, f"kmer{k}synteny"

def read_files(setname, synteny=None):
    """Read previously-calculated homology/synteny files and file frame."""
    set_path = Path(setname)
    files_frame_path = set_path/ f"{setname}{FILES_ENDING}"
    try:
        file_frame = pd.read_csv(files_frame_path, index_col=0, sep="\t")
    except:
        logger.error(f"Unable to read files frome from {files_frame_path}")
        sys.exit(1)
    if synteny is None:
        ending = HOMOLOGY_ENDING
        file_type = "homology"
    else:
        ending = f"-{synteny}{SYNTENY_ENDING}"
        file_type = "synteny"
    paths = [p for p in set_path.glob("*" + ending)]
    stems = [p.name[:-len(ending)] for p in paths]
    if len(stems) != len(file_frame):
        logger.error(f"Number of {file_type} files ({len(stems)})is not the same as length of file frame({len(file_frame)}).")
        sys.exit(1)
    frame_dict = {}
    for i, path in enumerate(paths):
        logger.debug(f"Reading homology file {path}.")
        frame_dict[stems[i]] = pd.read_csv(path, index_col=0, sep="\t")
    return file_frame, frame_dict


@cli.command()
@click_loguru.init_logger()
@click.option("--identity", "-i", default=0.0, help="Minimum sequence ID (0-1). [default: lowest]")
@click.option(
    "--clust/--no-clust",
    "-c/-x",
    is_flag=True,
    default=True,
    help="Do cluster calc.",
    show_default=True
)
@click.option("-s", "--shorten_source", default=False, is_flag=True, show_default=True, help="Remove invariant dotpaths in source IDs.")
@click.argument("setname")
@click.argument("gff_faa_path_list", nargs=-1)
def annotate_homology(identity, clust, shorten_source, setname, gff_faa_path_list):
    """Marshal homology and sequence information.

    Corresponding GFF and FASTA files must have a corresponding prefix to their
    file names, but theu may occur in any order in the list.  Paths to files
    need not be the same.  Files must be uncompressed. FASTA files must be
    protein files with extension ".faa".  GFF files must have extension ".gff3".

    IDs must correspond between GFF and FASTA files and must be unique across
    the entire set.
    """
    if not len(gff_faa_path_list):
        logger.error("No files in list, exiting.")
        sys.exit(0)
    # find matching pairs of gff files
    file_dict = {}
    gff_stems =  [str(Path(n).stem) for n in gff_faa_path_list if n.find(GFF_EXT) > -1]
    gff_stems.sort(key=lambda x:len(x))
    fasta_stems = [ str(Path(n).stem) for n in gff_faa_path_list if n.find(FASTA_EXT) > -1]
    fasta_stems.sort(key=lambda x:len(x))
    if len(gff_stems) != len(fasta_stems):
        logger.error(f"Differing number of {FASTA_EXT} ({len(fasta_stems)}) and {GFF_EXT} files ({len(gff_stems)}).")
        sys.exit(1)
    for fasta in fasta_stems:
        prefix_len = max([len(os.path.commonprefix([fasta, gff])) for gff in gff_stems])
        match_gff_idx = [i for i, gff in enumerate(gff_stems)
                     if len(os.path.commonprefix([fasta, gff])) == prefix_len][0]
        match_gff = gff_stems.pop(match_gff_idx)
        fasta_path = [Path(p) for p in gff_faa_path_list if p.endswith(fasta + "." + FASTA_EXT)][0]
        gff_path = [Path(p) for p in gff_faa_path_list if p.endswith(match_gff + "." +GFF_EXT)][0]
        stem = os.path.commonprefix([fasta, match_gff])
        file_dict[stem] = {GFF_EXT: gff_path,
                                                               FASTA_EXT: fasta_path}
    frame_dict = {}
    set_path = Path(setname)
    set_path.mkdir(parents=True, exist_ok=True)
    # TODO-all combinations of file sets
    fasta_records = []
    for stem in file_dict.keys():
        logger.debug(f"Reading GFF file {file_dict[stem][GFF_EXT]}.")
        annotation = gffpd.read_gff3(file_dict[stem][GFF_EXT])
        mRNAs = annotation.filter_feature_of_type(['mRNA']).attributes_to_columns()
        mRNAs.drop(mRNAs.columns.drop(['seq_id', 'start','strand', 'ID']), axis=1, inplace=True) # drop non-essential columns
        if shorten_source:
            # drop identical sub-fields in seq_id to keep them visually short (for development)
            split_sources = mRNAs["seq_id"].str.split(".", expand=True)
            split_sources = split_sources.drop([i for i in split_sources.columns if len(set(split_sources[i])) == 1], axis=1)
            sources = split_sources.agg(".".join, axis=1)
            mRNAs["seq_id"] = sources
        # TODO-calculated subfragments from repeated
        # TODO-sort GFFs in order of longest fragments
        # TODO-add gene order
        file_dict[stem]["fragments"] = len(set(mRNAs["seq_id"]))
        logger.debug(f"Reading FASTA file {file_dict[stem][FASTA_EXT]}.")
        fasta_dict = SeqIO.to_dict(SeqIO.parse(file_dict[stem][FASTA_EXT], "fasta"))
        # TODO-filter out crap and calculate ambiguous
        file_dict[stem]["n_seqs"] = len(fasta_dict)
        file_dict[stem]["residues"] = sum([len(fasta_dict[k].seq) for k in fasta_dict.keys()])
        mRNAs = mRNAs[mRNAs["ID"].isin(fasta_dict.keys())]
        mRNAs["protein_len"] = mRNAs["ID"].map(lambda k: len(fasta_dict[k].seq))
        frame_dict[stem]= mRNAs.set_index('ID')
        del annotation
        for key in fasta_dict.keys():
            fasta_records.append(fasta_dict[key])
    file_frame = pd.DataFrame.from_dict(file_dict).transpose()
    file_frame = file_frame.sort_values(by=["n_seqs"])
    file_frame["n"] = range(len(file_frame))
    file_frame["stem"] = file_frame.index
    file_frame = file_frame.set_index('n')
    logger.debug("Writing files frame.")
    file_frame.to_csv(set_path/ f"{setname}{FILES_ENDING}", sep="\t")
    del file_dict
    set_keys = list(file_frame["stem"])
    concatenated_fasta_name = f"{setname}.faa"
    if clust:
        logger.debug(f"Writing concatenated FASTA file {concatenated_fasta_name}.")
        with (set_path / concatenated_fasta_name).open("w") as concat_fh:
            SeqIO.write(fasta_records, concat_fh, "fasta")
        logger.debug("Doing cluster calculation.")
        cwd = Path.cwd()
        os.chdir(set_path)
        stats, graph, hist, any, all = usearch_cluster.callback(
            concatenated_fasta_name, identity, write_ids=True, delete=False)
        os.chdir(cwd)
        del stats, graph, hist, any, all
    del fasta_records
    cluster_frame = pd.read_csv(set_path/ (cluster_set_name(setname, identity) + "-ids.tsv"), sep="\t")
    cluster_dict = {}
    logger.debug("Getting ID-to-cluster info.")
    for index, row in cluster_frame.iterrows():
        if row['siz'] > 1:
            cluster_dict[row['id']] = row['cluster']
        else:
            cluster_dict[row['id']] = -1
    logger.debug("Mapping FASTA IDs to cluster.")
    for stem in set_keys:
        frame_dict[stem]["cluster"] = frame_dict[stem].index.map(cluster_dict)
    logger.debug(f"Writing homology frames.")
    for stem in set_keys:
        frame_dict[stem].to_csv(set_path / f"{stem}{HOMOLOGY_ENDING}", sep="\t")


@cli.command()
@click_loguru.init_logger()
@click.option("-k", default=6, help="Synteny block length.", show_default=True)
@click.option("-r", "--rmer", default=False, is_flag=True, show_default=True, help="Allow repeats in block.")
@click.argument("setname")
def kmer_synteny(k, rmer, setname):
    """Calculate k-mer syntenic blocks among sets of GFF/FASTA files.
    """
    set_path = Path(setname)
    files_frame, frame_dict = read_files(setname)
    set_keys = list(files_frame["stem"])
    logger.debug(f"Calculating k-mer of length {k} synteny blocks.")
    merge_frame_columns = ["hash", "source"]
    merge_frame = pd.DataFrame(columns=merge_frame_columns)
    for stem in set_keys:
        frame = frame_dict[stem]
        if rmer:
            map_func = None
            synteny_func_name = "None"
        else:
            map_func, synteny_func_name = kmer_synteny_block_func(k, frame)
        frame["footprint"] = 0
        frame["hashdir"] = 0
        frame[synteny_func_name] = 0
        frame_len = frame.shape[0]
        for idx in range(frame_len):
            footprint, hashdir, block = map_func()
            frame.iloc[idx, FOOTPRINT_POS] = footprint
            frame.iloc[idx, HASHDIR_POS] = hashdir
            frame.iloc[idx, SYNTENY_POS] = block
    #TODO:E values
        hash_series = frame.loc[:, synteny_func_name]
        assigned_hashes = hash_series[hash_series != 0]
        del hash_series
        n_assigned = len(assigned_hashes)
        logger.info(f"{stem} has {frame_len} proteins, {n_assigned} of which have {synteny_func_name} hashes,")
        hash_counts = assigned_hashes.value_counts()
        assigned_hash_frame  = pd.DataFrame(columns=merge_frame_columns)
        assigned_hash_frame["hash"] = assigned_hashes.unique()
        assigned_hash_frame["source"] = stem
        n_non_unique = n_assigned - len(assigned_hash_frame)
        percent_non_unique = n_non_unique/n_assigned * 100.
        logger.info(f"  of which {n_non_unique} ({percent_non_unique})% are non-unique.")
        merge_frame.append(assigned_hash_frame)
        del assigned_hash_frame
        # create self_count column in frame
        frame["self_count"] = 0
        for idx, row in frame[frame[synteny_func_name] != 0].iterrows():
            frame.loc[idx, "self_count"] = hash_counts.loc[row[synteny_func_name]]
        del hash_counts
    logger.debug(f"Calculating overlap of {len(merge_frame)} hash terms.")
    hash_counts = merge_frame["hash"].value_counts()
    merged_hash_frame = pd.DataFrame(index=merge_frame["hash"].unique(), columns=["count"])
    for idx, row in merged_hash_frame.iterrows():
        merged_hash_frame.loc[idx, "count"] = hash_counts.loc[row[synteny_func_name]]
    print(f"Merged_hash_frame={merged_hash_frame}")
    merged_hash_frame = merged_hash_frame[merged_hash_frame["count"] > 1]
    print(f"after dropping non-matching hashes, len = {len(merged_hash_frame)}")
    print(f"merged hash counts={hash_counts}")
    for stem in set_keys:
        synteny_name =  f"{stem}-{synteny_func_name}{SYNTENY_ENDING}"
        logger.debug(f"Writing {synteny_func_name} synteny frame {synteny_name}.")
        frame_dict[stem].to_csv(set_path / synteny_name, sep="\t")

def dagchainer_id_to_int(id):
    """Accepts DAGchainer ids such as "cl1" and returns an integer."""
    if not id.startswith("cl"):
        raise ValueError(f"Invalid ID {id}.")
    id_val = id[2:]
    if not id_val.isnumeric():
        raise ValueError(f"Non-numeric ID value in {id}.")
    return int(id_val)

@cli.command()
@click_loguru.init_logger()
@click.argument("setname")
def dagchainer_synteny(setname):
    """Read DAGchainer synteny blocks into homology frames.

    IDs must correspond between DAGchainer files and homology blocks.
    Currently does not calculate DAGchainer synteny.
    """
    synteny_func_name = "dagchainer"
    set_path = Path(setname)
    logger.debug(f"Reading {synteny_func_name} synteny file.")
    synteny_frame = pd.read_csv(set_path/ synteny_func_name / "clusters.tsv", sep="\t")
    synteny_frame["synteny_id"] = synteny_frame["cluster"].map(lambda id: dagchainer_id_to_int(id))
    synteny_frame = synteny_frame.drop(["cluster"], axis=1)
    cluster_counts = synteny_frame["synteny_id"].value_counts()
    synteny_frame["synteny_count"] = synteny_frame["synteny_id"].map(cluster_counts)
    synteny_frame = synteny_frame.sort_values(by=["synteny_count", "synteny_id"])
    synteny_frame = synteny_frame.set_index(["id"])
    files_frame, frame_dict = read_files(setname)
    set_keys = list(files_frame["stem"])
    def id_to_synteny_property(id, column):
        try:
            return synteny_frame.loc[id, column]
        except KeyError:
            return 0
    for stem in set_keys:
        homology_frame = frame_dict[stem]
        homology_frame["synteny_id"] = homology_frame.index.map(lambda x: id_to_synteny_property(x, "synteny_id"))
        homology_frame["synteny_count"] = homology_frame.index.map(lambda x: id_to_synteny_property(x, "synteny_count"))
        synteny_name =  f"{stem}-{synteny_func_name}{SYNTENY_ENDING}"
        logger.debug(f"Writing {synteny_func_name} synteny frame {synteny_name}.")
        homology_frame.to_csv(set_path / synteny_name, sep="\t")


@cli.command()
@click_loguru.init_logger()
@click.argument("setname")
@click.argument("synteny_type")
@click.argument("prefs", nargs=-1)
def proxy_genes(setname, synteny_type, prefs):
    """Calculate a set of proxy genes from synteny files.

    prefs is an optional list of genome stems in order of preference in the proxy calc.
    """
    set_path = Path(setname)
    files_frame, frame_dict = read_files(setname, synteny=synteny_type)
    set_keys = list(files_frame["stem"])
    if prefs != ():
        for stem in prefs:
            if stem not in set_keys:
                logger.error(f'Preference {stem} not in {set_keys}')
                sys.exit(1)
        logger.debug(f"Genome preference order in proxy selection is {prefs}")
    else:
        logger.debug("No preference list was set.")
    proxy_frame = None
    for stem in set_keys:
        frame_dict[stem]["stem"] = stem
        if proxy_frame is None:
            proxy_frame = frame_dict[stem]
        else:
            proxy_frame.append(frame_dict[stem])
    del files_frame
    proxy_frame = proxy_frame.sort_values(by=["cluster", "synteny_count", "synteny_id"])
    proxy_filename = f"{setname}-{synteny_type}{PROXY_ENDING}"
    logger.debug(f"Writing proxy file {proxy_filename}.")
    proxy_frame.to_csv(set_path/proxy_filename, sep="\t")

