# -*- coding: utf-8 -*-
"""External merges for synteny operations."""

# standard library imports
import array

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import networkx as nx
import numpy as np
import pandas as pd
from itertools import combinations
from loguru import logger


class AmbiguousMerger(object):

    """Counts instance of merges."""

    def __init__(
        self,
        graph_path=Path("synteny.gml"),
        count_key="value",
        ordinal_key="count",
        ambig_key="ambig",
        start_base=0,
        k=None,
        graph=None,
    ):
        """Create arrays as instance attributes."""
        self.graph_path = graph_path
        self.count_key = count_key
        self.ordinal_key = ordinal_key
        self.start_base = start_base
        self.ambig_key = ambig_key
        self.values = array.array("L")
        self.counts = array.array("h")
        self.ambig = array.array("h")
        if graph is None:
            self.graph = nx.Graph()
        else:
            self.graph = graph
        self.k = k
        self.n_edges = 0
        # self.print_counter = 0
        # self.print_increment = 1000

    def _unpack_payloads(self, vec):
        """Unpack TSV ints in payload"""
        wheres = np.where(~vec.mask)[0]
        values = np.array(
            [[int(i) for i in s.split("\t")] for s in vec.compressed()]
        ).transpose()
        return wheres, values

    def merge_func(self, value, count, payload_vec):
        """Return list of merged values."""
        self.values.append(value)
        self.counts.append(count)
        wheres, arr = self._unpack_payloads(payload_vec)
        max_ambig = arr[0].max()
        self.ambig.append(max_ambig)
        # if (self.print_counter%self.print_increment == 0):
        #    print(f"{value}: {max_ambig} {wheres} {arr}")
        # self.print_counter += 1
        if max_ambig == 1:
            self._adjacency_to_graph([f"{i}" for i in arr[1]], value)

    def _adjacency_to_graph(self, nodes, edgename):
        """Turn adjacency data into GML graph file."""
        self.graph.add_nodes_from(nodes)
        edges = list(combinations(nodes, 2))
        self.n_edges += len(edges)
        self.graph.add_edges_from(edges, label=f"{edgename}")

    def results(self):
        merge_frame = pd.DataFrame(
            {self.count_key: self.counts, self.ambig_key: self.ambig},
            index=self.values,
            dtype=pd.UInt32Dtype(),
        )
        merge_frame.sort_values(
            by=[self.ambig_key, self.count_key], inplace=True
        )
        unambig_frame = merge_frame[merge_frame[self.ambig_key] == 1].copy()
        n_unambig = len(unambig_frame)
        unambig_frame[self.ordinal_key] = (
            pd.array(
                range(self.start_base, self.start_base + n_unambig),
                dtype=pd.UInt32Dtype(),
            )
            * self.k
        )
        # merge_frame.to_csv('merge_frame.tsv', sep="\t")
        del unambig_frame[self.ambig_key]
        ambig_frame = merge_frame[merge_frame[self.ambig_key] > 1].copy()
        ambig_frame = ambig_frame.rename(
            columns={self.count_key: self.count_key + ".ambig"}
        )
        ambig_frame[self.ordinal_key + ".ambig"] = (
            pd.array(
                range(
                    self.start_base + n_unambig + 1,
                    self.start_base + len(ambig_frame) + n_unambig + 1,
                ),
                dtype=pd.UInt32Dtype(),
            )
            * self.k
        )
        # ambig_frame.to_csv('ambig_frame.tsv', sep="\t")
        del merge_frame
        component_sizes = [len(g) for g in nx.connected_components(self.graph)]
        # logger.info(f"Synteny graph has {self.graph.number_of_nodes()} nodes" +
        #            f" and {self.graph.number_of_edges()}/{self.n_edges}  edges")
        # logger.info(f"There are {len(component_sizes)} connected components of sizes {component_sizes}")
        nx.write_gml(self.graph, self.graph_path)
        return unambig_frame, ambig_frame, self.graph
