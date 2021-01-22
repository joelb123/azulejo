#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# stdlib imports
import sys
from pathlib import Path

# third-party imports
import click
import pandas as pd


@click.command()
@click.argument("synfile", type=click.Path(readable=True))
@click.argument("homfile", type=click.Path(readable=True))
def join_homology_synteny(synfile, homfile):
    "From a 2-column tsv INFILE, output a frequency-sorted file"
    synpath = Path(synfile)
    outpath = synpath.parent / (
        synpath.name[: -len(synpath.suffix)] + ".hom.tsv"
    )
    syn = pd.read_csv(synpath, sep="\t")
    hom = pd.read_csv(homfile, sep="\t")
    hom = hom.set_index("members")
    syn["hom_cluster"] = syn["members"].map(hom["cluster_id"])
    syn.to_csv(outpath, sep="\t")
    n_clusts = 0
    n_consistent = 0
    for cluster_id, subframe in syn.groupby(by=["cluster_id"]):
        n_clusts += 1
        hom_counts = subframe["cluster_id"].value_counts()
        if len(hom_counts) == 1:
            n_consistent += 1
    n_inconsistent = n_clusts - n_consistent
    print(f"   Inconsistent:\t{n_inconsistent}")


if __name__ == "__main__":
    join_homology_synteny()
