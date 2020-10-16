#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# stdlib imports
from pathlib import Path

# third-party imports
import click
import pandas as pd


@click.command()
@click.argument("infile", type=click.Path(readable=True))
@click.argument("outfile", type=click.Path(writable=True))
def rename_renumber_clusters(infile, outfile):
    "From a 2-column tsv INFILE, output a frequency-sorted file"
    outpath = Path(outfile)
    histpath = outpath.parent / (
        outpath.name[: -len(outpath.suffix)] + "_hist.tsv"
    )
    inclusts = pd.read_csv(infile, sep="\t", names=["cluster_id", "members"])
    cluster_counts = inclusts["cluster_id"].value_counts()
    cluster_map = pd.Series(
        range(len(cluster_counts)), index=cluster_counts.index
    )
    cluster_ids = inclusts["cluster_id"].map(cluster_map)
    cluster_sizes = inclusts["cluster_id"].map(cluster_counts)
    outclusts = pd.DataFrame(
        {
            "cluster_id": cluster_ids,
            "size": cluster_sizes,
            "members": inclusts["members"],
        }
    )
    outclusts.sort_values(by=["cluster_id"], inplace=True)
    outclusts.drop(
        outclusts[outclusts["size"] < 2].index, axis=0, inplace=True
    )
    n_clusts = len(outclusts)
    outclusts.index = range(n_clusts)
    print(f"Writing {n_clusts} clusters to {outpath}, histogram to {histpath}")
    outclusts.to_csv(outpath, sep="\t")
    cluster_hist = pd.DataFrame(cluster_counts.value_counts()).sort_index()
    cluster_hist = cluster_hist.rename(columns={"cluster_id": "n"})
    cluster_hist["cluster_pct"] = (
        cluster_hist["n"] * cluster_hist.index * 100.0 / n_clusts
    )
    cluster_hist.to_csv(histpath, sep="\t", float_format="%5.2f")


if __name__ == "__main__":
    rename_renumber_clusters()
