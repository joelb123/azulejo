#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# stdlib imports
from pathlib import Path

# third-party imports
import click
import networkx as nx
import pandas as pd


@click.command()
@click.argument("infile", type=click.Path(readable=True))
@click.argument("outfile", type=click.Path(writable=True))
def pairs_to_adjacency(infile, outfile):
    "From a 2-column tsv INFILE, output a frequency-sorted file"
    outpath = Path(outfile)
    histpath = outpath.parent / (
        outpath.name[: -len(outpath.suffix)] + "_hist.tsv"
    )
    G = nx.read_edgelist(infile)
    c = sorted(nx.connected_components(G), key=len, reverse=True)
    fh = outpath.open("w")
    fh.write("idx\tcluster_id\tsize\tmembers\n")
    line_no = 0
    count_list = []
    for i, comp in enumerate(c):
        size = len(comp)
        count_list.append(size)
        for node in comp:
            fh.write(f"{line_no}\t{i}\t{size}\t{node}\n")
            line_no += 1
    fh.close()
    print(f"{line_no} items in adjacency written to {outpath}")
    del G, c
    n_clusts = len(count_list)
    n_items = sum(count_list)
    cluster_counts = pd.DataFrame({"size": count_list})
    cluster_hist = (
        pd.DataFrame(cluster_counts.value_counts()).sort_index().reset_index()
    )
    cluster_hist = cluster_hist.set_index("size")
    cluster_hist = cluster_hist.rename(columns={0: "n"})
    cluster_hist["item_pct"] = (
        cluster_hist["n"] * cluster_hist.index * 100.0 / n_items
    )
    cluster_hist.to_csv(histpath, sep="\t", float_format="%5.2f")
    cluster_hist["cluster_pct"] = cluster_hist["n"] * 100.0 / n_clusts
    cluster_hist.to_csv(histpath, sep="\t", float_format="%5.2f")


if __name__ == "__main__":
    pairs_to_adjacency()
