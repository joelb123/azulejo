# -*- coding: utf-8 -*-
"""Synteny (genome order) operations."""
# standard library imports
import statistics
import sys

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import numpy as np
import pandas as pd

# first-party imports
from loguru import logger


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
                f"subcluster {subcluster} is not unique w.r.t. proteome"
                f" {list(stems)[best_choice]}."
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
            f"Dropping {len(self.drop_ids)} ({drop_pct:0.1f}%) of"
            f" {len(self.frame)} genes."
        )
        return self.frame.drop(self.drop_ids)

    def selection_stats(self):
        """Return selection stats."""
        return (
            self.cluster_count,
            self.first_choice_unavailable,
            self.first_choice_hits,
        )


def calculate_proxy_genes(setname, synteny_type, prefs):
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
        by=["cluster_size", "hom.cluster", "synteny_count", "synteny_id"]
    )
    proxy_filename = f"{setname}-{synteny_type}-{PROXY_ENDING}"
    logger.debug(f"Writing initial proxy file {proxy_filename}.")
    proxy_frame.to_csv(set_path / proxy_filename, sep="\t")
    proxy_frame["reason"] = ""
    logger.debug("Downselecting homology clusters.")
    downselector = ProxySelector(proxy_frame, prefs)
    for unused_cluster, homology_cluster in proxy_frame.groupby(
        by=["hom.cluster"]
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
        f"First-choice ({prefs[0]}) selections from {cluster_count} homology"
        " clusters:"
    )
    logger.info(
        f"   not in cluster: {first_choice_unavailable}"
        f" ({first_choice_unavailable_percent:.1f}%)"
    )
    logger.info(
        f"   chosen as proxy: {first_choice_hits}"
        f" ({first_choice_percent:.1f}%)"
    )
