# -*- coding: utf-8 -*-
"""Send and receive on-disk messages through names pips with file locking."""
# standard library imports
import array
import contextlib
import fcntl
from pathlib import Path

# third-party imports
import attr
import numpy as np
import numpy.ma as ma
import pandas as pd


@attr.s
class DataMailboxes:

    """Pass data to and from on-disk FIFOs."""

    n_boxes = attr.ib()
    mb_dir_path = attr.ib(default=Path("./mailboxes/"))
    file_extension = attr.ib(default=None)

    def write_headers(self, header):
        """Initialize the mailboxes, optionally writing header."""
        self.mb_dir_path.mkdir(parents=True, exist_ok=True)
        for i in range(self.n_boxes):
            mb_path = self.path_to_mailbox(i)
            with mb_path.open("w") as fh:
                fh.write(header)

    @contextlib.contextmanager
    def locked_open_for_write(self, box_no):
        """Open a mailbox with an advisory file lock."""
        mb_path = self.path_to_mailbox(box_no)
        with mb_path.open("a+") as fd:
            fcntl.flock(fd, fcntl.LOCK_EX)
            yield fd
            fcntl.flock(fd, fcntl.LOCK_UN)

    @contextlib.contextmanager
    def open_then_delete(self, box_no, delete=True):
        """Open a mailbox with an advisory file lock."""
        box_path = self.path_to_mailbox(box_no)
        with box_path.open("r") as fd:
            yield fd
            if delete:
                box_path.unlink()

    def path_to_mailbox(self, box_no):
        """Return a path to a mailbox file."""
        if self.file_extension is None:
            ext = ""
        else:
            ext = f".{self.file_extension}"
        return self.mb_dir_path / f"{box_no}{ext}"

    def delete(self):
        """Remove the mailbox directory. """
        file_list = list(self.mb_dir_path.glob("*"))
        for file in file_list:
            file.unlink()
        self.mb_dir_path.rmdir()


@attr.s
class ExternalMerge(object):

    """Merges integers from files."""

    file_path_func = attr.ib(default=None)
    n_merge = attr.ib(default=None)
    value_vec = None
    fh_list = []

    def init(self, header_value):
        """Read and check the header info."""
        self.fh_list = [
            self.file_path_func(i).open() for i in range(self.n_merge)
        ]
        assert [next(fh).rstrip() for fh in self.fh_list] == [
            header_value
        ] * self.n_merge
        self.value_vec = ma.masked_array(
            np.zeros(self.n_merge), mask=np.zeros(self.n_merge),
        ).astype(np.uint64)
        self._next_vals(np.arange(self.n_merge))

    def _next_vals(self, index_vec):
        """Return the value of the next iterated value."""
        for index in index_vec:
            try:
                self.value_vec[index] = int(next(self.fh_list[index]))
            except StopIteration:
                self.value_vec.mask[index] = 1

    def _close_all(self):
        """Close all filehandles."""
        for i in range(self.n_merge):
            self.fh_list[i].close()

    def merge(self, count_key=None, ordinal_key=None):
        """Return list of merged values."""
        hashes = array.array("L")
        counts = array.array("h")
        while (~self.value_vec.mask).sum() > 1:
            minimum = np.amin(self.value_vec)
            where_min = np.where(self.value_vec == minimum)[0]
            self._next_vals(where_min)
            count = len(where_min)
            if count > 1:
                hashes.append(minimum)
                counts.append(count)
        self._close_all()
        merge_frame = pd.DataFrame(
            pd.array(counts, dtype=pd.UInt32Dtype()),
            columns=[count_key],
            index=hashes,
        )
        merge_frame.sort_values(by=[count_key], inplace=True)
        merge_frame["tmp.base"] = pd.array(
            range(len(merge_frame)), dtype=pd.UInt32Dtype()
        )
        return merge_frame
