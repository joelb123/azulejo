# -*- coding: utf-8 -*-
"""Synteny hash operations."""


# third-party imports
import attr
import numpy as np
import pandas as pd

# module imports
from .common import DIRECTIONAL_CATEGORY
from .common import hash_array
from .common import logger

# helper functions


def _cum_val_count(arr):
    """Return an array of cumulative counts of values."""
    counts = {}
    out_arr = pd.array([pd.NA] * len(arr), dtype=pd.UInt32Dtype(),)
    for i, val in enumerate(arr):
        if pd.isnull(val):
            continue
        elif val in counts:
            counts[val] += 1
        else:
            counts[val] = 1
        out_arr[i] = counts[val]
    return out_arr


def _true_positions_and_runs(bool_vec):
    """Return arrays of positions and lengths of runs of True."""
    uneq_idxs = np.append(
        np.where(bool_vec[1:] != bool_vec[:-1]), bool_vec.size - 1
    )
    runlengths = np.diff(np.append(-1, uneq_idxs))
    positions = np.cumsum(np.append(0, runlengths))[:-1]
    true_idxs = np.where(bool_vec[positions])
    return positions[true_idxs], runlengths[true_idxs]


def _fill_na_with_last_valid(ser, flip=False):
    """Input a series with NA values, returns a series with those values filled."""
    lv_arr = pd.array([pd.NA] * len(ser), dtype=pd.UInt32Dtype(),)
    if not (ser.isnull().all() or ser.notna().all()):
        null_vec = ser.isnull().to_numpy()
        val_vec = ser.to_numpy()
        if flip:
            null_vec = np.flip(null_vec)
            val_vec = np.flip(val_vec)
        first_null_pos, null_runs = _true_positions_and_runs(null_vec)
        fill_vals = np.append(pd.NA, val_vec)[first_null_pos]
        for i, pos in enumerate(first_null_pos):
            for j in range(null_runs[i]):
                lv_arr[pos + j] = fill_vals[i]
        if flip:
            lv_arr = np.flip(lv_arr)
        lv_ser = pd.Series(lv_arr, index=ser.index)
        return lv_ser


def _cum_val_cnt_where_ser2_is_na(ser1, ser2, flip=False):
    """Return the cumulative value count of ser1 in regions where ser2 is NA."""
    if len(ser1) != len(ser2):
        logger.warning(f"Lengths of ser1 and ser2 differ at {ser1}")
    vc_arr = pd.array([pd.NA] * len(ser1), dtype=pd.UInt32Dtype(),)
    if not (ser2.isnull().all() or ser2.notna().all()):
        null_vec = ser2.isnull().to_numpy()
        val_vec = ser1.to_numpy()
        if flip:
            null_vec = np.flip(null_vec)
            val_vec = np.flip(val_vec)
        null_pos, null_runs = _true_positions_and_runs(null_vec)
        null_len = len(null_pos)
        for i in range(null_len):
            vc_arr[
                null_pos[i] : (null_pos[i] + null_runs[i])
            ] = _cum_val_count(
                val_vec[null_pos[i] : (null_pos[i] + null_runs[i])]
            )
        if flip:
            vc_arr = np.flip(vc_arr)
    vc_ser = pd.Series(vc_arr, index=ser2.index)
    return vc_ser


def calculate_disambig_hashes(df):
    """Calculate disambiguation frame (per-fragment)."""
    hash2_fr = df[["syn.anchor.id", "tmp.ambig.id"]].copy()
    hash2_fr = hash2_fr.rename(columns={"syn.anchor.id": "tmp.anchor.id"})
    hash2_fr["tmp.upstr_anchor"] = _fill_na_with_last_valid(
        df["syn.anchor.id"]
    )
    hash2_fr["tmp.downstr_anchor"] = _fill_na_with_last_valid(
        df["syn.anchor.id"], flip=True
    )
    hash2_fr["tmp.upstr_occur"] = _cum_val_cnt_where_ser2_is_na(
        df["tmp.ambig.id"], df["syn.anchor.id"]
    )
    hash2_fr["tmp.downstr_occur"] = _cum_val_cnt_where_ser2_is_na(
        df["tmp.ambig.id"], df["syn.anchor.id"], flip=True
    )
    hash2_fr["tmp.i"] = range(len(hash2_fr))
    upstream_hash = pd.array([pd.NA] * len(hash2_fr), dtype=pd.UInt32Dtype())
    downstream_hash = pd.array([pd.NA] * len(hash2_fr), dtype=pd.UInt32Dtype())
    hash2_fr["tmp.disambig.up"] = pd.NA
    hash2_fr["tmp.disambig.down"] = pd.NA
    for unused_id, row in hash2_fr.iterrows():
        row_no = row["tmp.i"]
        ambig_base = row["tmp.ambig.id"]
        upstream_unambig = row["tmp.upstr_anchor"]
        downstream_unambig = row["tmp.downstr_anchor"]
        occur_upstream = row["tmp.upstr_occur"]
        occur_downstream = row["tmp.downstr_occur"]
        if pd.notna(ambig_base):
            if pd.notna(upstream_unambig):
                if not pd.notna(occur_upstream):
                    logger.warning(
                        f"Something is wrong upstream of base {ambig_base}"
                    )
                upstream_hash[row_no] = hash_array(
                    np.array([upstream_unambig, ambig_base, occur_upstream])
                )
            if pd.notna(downstream_unambig):
                if not pd.notna(occur_downstream):
                    logger.warning(
                        f"Something is wrong downstream of base {ambig_base}"
                    )
                downstream_hash[row_no] = hash_array(
                    np.array(
                        [ambig_base, downstream_unambig, occur_downstream]
                    )
                )
    hash2_fr["tmp.disambig.up"] = upstream_hash
    hash2_fr["tmp.disambig.down"] = downstream_hash
    return hash2_fr[["tmp.disambig.up", "tmp.disambig.down"]]


@attr.s
class SyntenyBlockHasher(object):
    """Synteny-block hashes via reversible-peatmer method."""

    k = attr.ib(default=3)
    peatmer = attr.ib(default=True)
    prefix = attr.ib(default="syn")

    def hash_name(self, no_prefix=False):
        """Return the string name of the hash function."""
        if no_prefix:
            prefix_str = ""
        else:
            prefix_str = self.prefix + "."
        if self.peatmer:
            return f"{prefix_str}hash.peatmer{self.k}"
        return f"{prefix_str}hash.kmer{self.k}"

    def shingle(self, cluster_series, direction, hash_val):
        """Return a vector of anchor ID's. """
        vec = cluster_series.to_numpy().astype(int)
        steps = np.insert((vec[1:] != vec[:-1]).astype(int), 0, 0).cumsum()
        if max(steps) != self.k - 1:
            logger.warning(f"Inconsistency in shingling hash {hash_val}")
            logger.warning(f"input homology string={vec}")
            steps[np.where(steps > self.k - 1)] = self.k - 1
        if direction == "-":
            steps = self.k - 1 - steps
        return steps

    def calculate(self, cluster_series):
        """Return an array of synteny block hashes data."""
        # Maybe the best code I've ever written--JB
        vec = cluster_series.to_numpy().astype(int)
        if self.peatmer:
            uneq_idxs = np.append(np.where(vec[1:] != vec[:-1]), vec.size - 1)
            runlengths = np.diff(np.append(-1, uneq_idxs))
            positions = np.cumsum(np.append(0, runlengths))[:-1]
            n_mers = len(positions) - self.k + 1
            footprints = pd.array(
                [runlengths[i : i + self.k].sum() for i in range(n_mers)],
                dtype=pd.UInt32Dtype(),
            )
        else:
            n_elements = len(cluster_series)
            n_mers = n_elements - self.k + 1
            positions = np.arange(n_elements)
            footprints = pd.array([self.k] * n_mers, dtype=pd.UInt32Dtype())
        if n_mers < 1:
            return None
        # Calculate k-mers over indirect index
        kmer_mat = np.array(
            [vec[positions[i : i + self.k]] for i in range(n_mers)]
        )
        fwd_rev_hashes = np.array(
            [
                np.apply_along_axis(hash_array, 1, kmer_mat),
                np.apply_along_axis(hash_array, 1, np.flip(kmer_mat, axis=1)),
            ]
        )
        plus_minus = np.array([["+"] * n_mers, ["-"] * n_mers])
        directions = np.take_along_axis(
            plus_minus,
            np.expand_dims(fwd_rev_hashes.argmin(axis=0), axis=0),
            axis=0,
        )[0]
        return pd.DataFrame(
            [
                pd.Categorical(directions, dtype=DIRECTIONAL_CATEGORY),
                footprints,
                pd.array(
                    np.amin(fwd_rev_hashes, axis=0), dtype=pd.UInt32Dtype()
                ),
            ],
            columns=[
                "syn.hash.direction",
                "syn.hash.footprint",
                self.hash_name(),
            ],
            index=cluster_series.index[positions[:n_mers]],
        )
