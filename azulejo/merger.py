# -*- coding: utf-8 -*-
"""External merges for synteny operations."""

# standard library imports
import array

# from os.path import commonprefix as prefix
from pathlib import Path

# third-party imports
import numpy as np
import pandas as pd

# helper_functions
def _unpack_payloads(vec):
    """Unpack TSV ints in payload."""
    values = np.array(
        [[int(i) for i in s.split("\t")] for s in vec.compressed()]
    ).transpose()
    return values


class AmbiguousMerger(object):
    """Counts instance of merges."""

    def __init__(
        self,
        count_key="value",
        ordinal_key="count",
        ambig_key="ambig",
        ambig_count_key=None,
        start_base=0,
        ambig_ordinal_key="count.ambig",
        alt_hash=False,
    ):
        """Create arrays as instance attributes."""
        self.count_key = count_key
        self.ordinal_key = ordinal_key
        self.start_base = start_base
        self.ambig_key = ambig_key
        self.ambig_ordinal_key = ambig_ordinal_key
        self.values = array.array("L")
        self.counts = array.array("h")
        self.ambig = array.array("h")
        self.alt_hash = alt_hash
        self.ambig_count_key = ambig_count_key
        if self.alt_hash:
            self.alt_hash_dict = {}
            self.count_dict = {}
            self.ambig_dict = {}
            self.id_no = self.start_base

    def merge_func(self, value, count, payload_vec):
        """Return list of merged values."""
        self.values.append(value)
        self.counts.append(count)
        arr = _unpack_payloads(payload_vec)
        max_ambig = arr[0].max()
        self.ambig.append(max_ambig)
        if self.alt_hash:
            # if this is true, vec[1] contains a list
            # of alt hashes, with 0 as a fill value to be ignoreed
            self.alt_hash_dict[value] = list({i for i in arr[1] if i != 0})
            self.count_dict[value] = count
            self.ambig_dict[value] = max_ambig

    def results(self):
        """Calculate list of merges."""
        drop_list = []
        if self.alt_hash:
            for hash_val in self.alt_hash_dict:
                related_hashes = [hash_val] + [
                    alt
                    for alt in self.alt_hash_dict[hash_val]
                    if alt in self.count_dict
                ]
                if len(related_hashes) == 1:
                    continue
                related_hashes.sort()  # take the first in numberical order if all else is equal
                non_ambig_hashes = [
                    h for h in related_hashes if self.ambig_dict[h] == 1
                ]
                max_count_idx = np.argmax(
                    [self.count_dict[h] for h in non_ambig_hashes]
                )
                best_hash = non_ambig_hashes[max_count_idx]
                if best_hash != hash_val:
                    drop_list.append(hash_val)
            del self.alt_hash_dict, self.count_dict, self.ambig_dict
        merge_frame = pd.DataFrame(
            {self.count_key: self.counts, self.ambig_key: self.ambig},
            index=self.values,
            dtype=pd.UInt32Dtype(),
        )
        merge_frame.drop(drop_list, inplace=True)
        merge_frame.sort_values(
            by=[self.ambig_key, self.count_key], inplace=True
        )
        unambig_frame = merge_frame[merge_frame[self.ambig_key] == 1].copy()
        n_unambig = len(unambig_frame)
        unambig_frame[self.ordinal_key] = pd.array(
            range(self.start_base, self.start_base + n_unambig),
            dtype=pd.UInt32Dtype(),
        )
        del unambig_frame[self.ambig_key]
        ambig_frame = merge_frame[merge_frame[self.ambig_key] > 1].copy()
        del merge_frame
        if self.ambig_count_key is None:
            # Don't pass counts along
            del ambig_frame[self.count_key]
        else:
            ambig_frame = ambig_frame.rename(
                columns={self.count_key: self.ambig_count_key}
            )
        del ambig_frame[self.ambig_key]
        ambig_frame[self.ambig_ordinal_key] = pd.array(
            range(
                self.start_base + n_unambig,
                self.start_base + len(ambig_frame) + n_unambig,
            ),
            dtype=pd.UInt32Dtype(),
        )
        return unambig_frame, ambig_frame
