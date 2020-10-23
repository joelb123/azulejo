#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# stdlib imports
from pathlib import Path

# third-party imports
import click
import xxhash
import networkx as nx
import numpy as np
import pandas as pd


def hash_array(kmer):
    """Return a hash of a numpy array."""
    return xxhash.xxh32_intdigest(kmer.tobytes())


@click.command()
@click.argument("infile", type=click.Path(readable=True))
@click.argument("outfile", type=click.Path(writable=True))
def pairs_to_adjacency(infile, outfile):
    "From a 2-column tsv INFILE, output a frequency-sorted file"
    outpath = Path(outfile)
    summarypath = outpath.parent / (
        outpath.name[: -len(outpath.suffix)] + "_summary.tsv"
    )
    histpath = outpath.parent / (
        outpath.name[: -len(outpath.suffix)] + "_hist.tsv"
    )
    G = nx.read_edgelist(infile)
    c = sorted(nx.connected_components(G), key=len, reverse=True)
    fh = outpath.open("w")
    fh.write("idx\tcluster_id\tsize\tmembers\n")
    n_items = 0
    count_list = []
    hash_list = []
    id_list = []
    for i, comp in enumerate(c):
        id_list.append(i)
        size = len(comp)
        count_list.append(size)
        hash_list.append(hash_array(np.array(sorted(comp))))
        for node in comp:
            fh.write(f"{n_items}\t{i}\t{size}\t{node}\n")
            n_items += 1
    fh.close()
    print(f"   Items:\t{n_items}")
    n_clusts = len(count_list)
    print(f"   Clusters:\t{n_clusts}")
    del G, c
    cluster_counts = pd.DataFrame({"size": count_list})
    largest_cluster = cluster_counts["size"].max()
    print(f"   Largest:\t{largest_cluster}")
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
    clusters = pd.DataFrame(
        {"anchor.id": id_list, "count": count_list, "hash": hash_list}
    )
    clusters.to_csv(summarypath, sep="\t")


if __name__ == "__main__":
    pairs_to_adjacency()
