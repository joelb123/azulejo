# -*- coding: utf-8 -*-
#
# standard library imports
#
import contextlib
import os
import sys
from collections import Counter, OrderedDict
from datetime import datetime
from itertools import chain, combinations

#
# third-party imports
#
import click
import dask.bag as db
from dask.diagnostics import ProgressBar
import gffpandas.gffpandas as gffpd
import networkx as nx
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data import IUPACData
from loguru import logger

#
# package imports
#
from . import cli
from . import click_loguru
from .common import *

# global constants
GFF_EXT = "gff3"
FASTA_EXT = "faa"

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
@click.argument("setname")
@click.argument("gff_faa_path_list", nargs=-1)
def synteny(identity, clust, setname, gff_faa_path_list):
    """Calculate syntenic blocks among sets of GFF/FASTA files.

    Filenames must correspond between GFF and FASTA files, but the may occur
    in any order in the list.  Files must be uncompressed.  Paths to files
    need not be the same. FASTA files must be protein files with extension
     ".faa".  GFF files must have extension ".gff3".

    IDs must correspond between GFF and FASTA files and must be unique across
    the entire set.
    """
    if not len(gff_faa_path_list):
        logger.error("No files in list, exiting.")
        sys.exit(0)
    filepaths = [Path(p) for p in gff_faa_path_list]
    filestems = set([f.stem for f in filepaths])
    file_dict = {}
    for stem in filestems:
        file_dict[stem] = {}
    for name in filestems:
        faa_paths = [p for p in filepaths if p.suffix.endswith(FASTA_EXT) and p.stem == name]
        gff_paths = [p for p in filepaths if p.suffix.endswith(GFF_EXT) and p.stem == name]
        if len(faa_paths) != 1 or len(gff_paths) != 1:
            logger.error(f"The wrong number of files with name {name} were found.")
            sys.exit(1)
        file_dict[name][FASTA_EXT] = faa_paths[0]
        file_dict[name][GFF_EXT] = gff_paths[0]
    frame_dict = {}
    for stem in filestems:
        logger.debug(f"Reading GFF file {file_dict[stem][GFF_EXT]}.")
        annotation = gffpd.read_gff3(file_dict[stem][GFF_EXT])
        mRNAs = annotation.filter_feature_of_type(['mRNA']).attributes_to_columns()
        mRNAs.drop(mRNAs.columns.drop(['seq_id', 'start','strand', 'ID']), axis=1, inplace=True)
        split_sources = mRNAs["seq_id"].str.split(".", expand=True)
        mRNAs["seq_id"] = split_sources.drop([i for i in split_sources.columns if len(set(split_sources[i])) == 1], axis=1).agg(".".join, axis=1)
        fasta_dict = SeqIO.to_dict(SeqIO.parse(file_dict[stem][FASTA_EXT], "fasta"))
        mRNAs = mRNAs[mRNAs["ID"].isin(fasta_dict.keys())]
        mRNAs["protein_len"] = mRNAs["ID"].map(lambda k: len(fasta_dict[k].seq))
        frame_dict[stem]= mRNAs.set_index('ID')
        del annotation
    #
    # call to usearch_cluster goes here
    #
    cluster_frame = pd.read_csv("all-nr-0000-ids.tsv", sep="\t")
    cluster_dict = {}
    for index, row in cluster_frame.iterrows():
        if row['siz'] > 1:
            cluster_dict[row['id']] = row['cluster']
        else:
            cluster_dict[row['id']] = None
    logger.debug("Mapping")
    for stem in filestems:
        frame_dict[stem]["cluster"] = frame_dict[stem].index.map(cluster_dict)
    #
    # calculate synteny blocks here
    #
    for stem in filestems:
        frame_dict[stem].to_csv(f"{stem}-synteny.tsv", sep="\t")
