# -*- coding: utf-8 -*-
"""Synteny hash operations."""


# third-party imports
import attr
import numpy as np
import pandas as pd
import xxhash
from loguru import logger

# module imports
from .common import DIRECTIONAL_CATEGORY
from .common import remove_tmp_columns


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

    def _hash_array(self, kmer):
        """Return a hash of an array."""
        return xxhash.xxh32_intdigest(kmer.tobytes())

    def shingle(self, cluster_series, base, direction, hash):
        """Return a vector of anchor ID's. """
        vec = cluster_series.to_numpy().astype(int)
        steps = np.insert((vec[1:] != vec[:-1]).astype(int), 0, 0).cumsum()
        try:
            assert max(steps) == self.k - 1
        except AssertionError:
            logger.error(
                f"Working around minor error in shingling hash {hash}, base"
                f" {base};"
            )
            logger.error(f"input homology string={vec}")
            logger.error(f"max index = {max(steps)}, should be {self.k-1}")
            steps[np.where(steps > self.k - 1)] = self.k - 1
        if direction == "+":
            return base + steps
        return base + self.k - 1 - steps

    def calculate_disambig_hashes(self, df):
        """Calculated Disambiguation frame (per-fragment)."""
        hash2_fr = df[["syn.anchor_id", "tmp.base.ambig"]].copy()
        hash2_fr = hash2_fr.rename(columns={"syn.anchor_id": "tmp.anchor_id"})
        hash2_fr["tmp.upstr_anchor"] = self.fill_na_with_last_valid(
            df["syn.anchor_id"]
        )
        hash2_fr["tmp.downstr_anchor"] = self.fill_na_with_last_valid(
            df["syn.anchor_id"], flip=True
        )
        hash2_fr["tmp.upstr_occur"] = self.cum_val_count_where_ser2_is_NA(
            df["tmp.base.ambig"], df["syn.anchor_id"]
        )
        hash2_fr["tmp.downstr_occur"] = self.cum_val_count_where_ser2_is_NA(
            df["tmp.base.ambig"], df["syn.anchor_id"], flip=True
        )
        hash2_fr["tmp.i"] = range(len(hash2_fr))
        upstream_hash = pd.array(
            [pd.NA] * len(hash2_fr), dtype=pd.UInt32Dtype()
        )
        downstream_hash = pd.array(
            [pd.NA] * len(hash2_fr), dtype=pd.UInt32Dtype()
        )
        hash2_fr["syn.disambig_upstr"] = pd.NA
        hash2_fr["syn.disambig_downstr"] = pd.NA
        for unused_id, row in hash2_fr.iterrows():
            row_no = row["tmp.i"]
            ambig_base = row["tmp.base.ambig"]
            upstream_unambig = row["tmp.upstr_anchor"]
            downstream_unambig = row["tmp.downstr_anchor"]
            occur_upstream = row["tmp.upstr_occur"]
            occur_downstream = row["tmp.downstr_occur"]
            if pd.notna(ambig_base):
                if pd.notna(upstream_unambig):
                    assert pd.notna(occur_upstream)
                    upstream_hash[row_no] = self._hash_array(
                        np.array(
                            [upstream_unambig, ambig_base, occur_upstream]
                        )
                    )
                if pd.notna(downstream_unambig):
                    assert pd.notna(occur_downstream)
                    downstream_hash[row_no] = self._hash_array(
                        np.array(
                            [ambig_base, downstream_unambig, occur_downstream]
                        )
                    )
        hash2_fr["syn.disambig_upstr"] = upstream_hash
        hash2_fr["syn.disambig_downstr"] = downstream_hash
        return remove_tmp_columns(hash2_fr)

    def fill_na_with_last_valid(self, ser, flip=False):
        """Input a series with NA values, returns a series with those values filled."""
        lv_arr = pd.array([pd.NA] * len(ser), dtype=pd.UInt32Dtype(),)
        if not (ser.isnull().all() or ser.notna().all()):
            null_vec = ser.isnull().to_numpy()
            val_vec = ser.to_numpy()
            if flip:
                null_vec = np.flip(null_vec)
                val_vec = np.flip(val_vec)
            first_null_pos, null_runs = self._true_positions_and_runs(null_vec)
            fill_vals = np.append(pd.NA, val_vec)[first_null_pos]
            for i, pos in enumerate(first_null_pos):
                for j in range(null_runs[i]):
                    lv_arr[pos + j] = fill_vals[i]
            if flip:
                lv_arr = np.flip(lv_arr)
        lv_ser = pd.Series(lv_arr, index=ser.index)
        return lv_ser

    def _true_positions_and_runs(self, bool_vec):
        """Return arrays of positions and lengths of runs of True."""
        runlengths, positions = self._run_lengths_and_positions(bool_vec)
        true_idxs = np.where(bool_vec[positions])
        return positions[true_idxs], runlengths[true_idxs]

    def cum_val_count_where_ser2_is_NA(self, ser1, ser2, flip=False):
        """Return the cumulative value count of ser1 in regions where ser2 is NA."""
        assert len(ser1) == len(ser2)
        vc_arr = pd.array([pd.NA] * len(ser1), dtype=pd.UInt32Dtype(),)
        if not (ser2.isnull().all() or ser2.notna().all()):
            null_vec = ser2.isnull().to_numpy()
            val_vec = ser1.to_numpy()
            if flip:
                null_vec = np.flip(null_vec)
                val_vec = np.flip(val_vec)
            null_pos, null_runs = self._true_positions_and_runs(null_vec)
            for i in range(len(null_pos)):
                vc_arr[
                    null_pos[i] : (null_pos[i] + null_runs[i])
                ] = self._cum_val_count(
                    val_vec[null_pos[i] : (null_pos[i] + null_runs[i])]
                )
            if flip:
                vc_arr = np.flip(vc_arr)
        vc_ser = pd.Series(vc_arr, index=ser2.index)
        return vc_ser

    def _cum_val_count(self, arr):
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

    def _run_lengths_and_positions(self, vec):
        """Compute vector of runlengths."""
        uneq_idxs = np.append(np.where(vec[1:] != vec[:-1]), len(vec) - 1)
        runlengths = np.diff(np.append(-1, uneq_idxs))
        positions = np.cumsum(np.append(0, runlengths))[:-1]
        return runlengths, positions

    def calculate(self, cluster_series):
        """Return an array of synteny block hashes data."""
        # Maybe the best code I've ever written--JB
        vec = cluster_series.to_numpy().astype(int)
        if self.peatmer:
            uneq_idxs = np.append(np.where(vec[1:] != vec[:-1]), len(vec) - 1)
            runlengths, positions = self._run_lengths_and_positions(vec)
            n_mers = len(positions) - self.k + 1
            footprints = pd.array(
                [runlengths[i : i + self.k].sum() for i in range(n_mers)],
                dtype=pd.UInt32Dtype(),
            )
        else:
            n_elements = len(cluster_series)
            n_mers = n_elements - self.k + 1
            positions = np.arange(n_elements)
            footprints = pd.array([self.k] * (n_mers), dtype=pd.UInt32Dtype())
        if n_mers < 1:
            return None
        # Calculate k-mers over indirect index
        kmer_mat = np.array(
            [vec[positions[i : i + self.k]] for i in range(n_mers)]
        )
        fwd_rev_hashes = np.array(
            [
                np.apply_along_axis(self._hash_array, 1, kmer_mat),
                np.apply_along_axis(
                    self._hash_array, 1, np.flip(kmer_mat, axis=1)
                ),
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
