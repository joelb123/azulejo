# -*- coding: utf-8 -*-
"""Send and receive on-disk messages through names pips with file locking."""
# standard library imports
import contextlib
import fcntl
import sys
from pathlib import Path

# third-party imports
import attr
import numpy as np
import numpy.ma as ma
from memory_tempfile import MemoryTempfile

# module imports
from .common import SCRATCH_DEV
from .common import append_slash
from .common import disk_usage_mb
from .common import free_mb
from .common import is_writable

# global constants
MAX_LINE_LEN = 144


# shared functions


class SpaceAwareTempDevices:
    """Allocate and collect stats on the fastest device that will hold your data."""

    class TempDirAllocator:
        """Reserve space and collect usage stats on a space-aware temp device."""

        def __init(self, space_needed_mb=None):
            """Select the fastest device with enough memory."""
            self.space_needed = space_needed_mb
            self.dev = fallback_dir
            if self.space_needed is not None:
                if space_needed_mb < FREE_MEM_MB:
                    memdev_list = MemoryTempfile(
                        filesystem_types=["tmpfs", "shm"]
                    ).get_usable_mem_tempdir_paths()
                    for memdev in memdev_list:
                        self.initial_space_free = free_mb(memdev)
                        if self.initial_space_free >= space_needed_mb:
                            self.dev = memdev
                            break
                if self.dev == fallback_dir and SCRATCH_DEV is not None:
                    self.initial_space_free = free_mb(SCRATCH_DEV)
                    if self.initial_space_free >= space_needed_mb:
                        self.dev = SCRATCH_DEV
            if self.dev == fallback_dir:
                self.initial_space_free = free_mb(self.dev)

    def __init__(
        self, scratch_dev=None, fallback_dev="/tmp/", addl_fs_types=None
    ):
        """Create a list of writable tmp devices and space allocations on them."""
        self.fallback_dev = fallback_dev
        self.scratch_dev = scratch_dev
        fs_types = ["tmpfs", "shm"]
        if addl_fs_types is not None:
            fs_types += addl_fs_types
        self.memdev_list = [
            append_slash(dev)
            for dev in MemoryTempfile(
                filesystem_types=fs_types
            ).get_usable_mem_tempdir_paths()
            if is_writable(dev)
        ]
        logger.debug(f"Memory devices are {self.memdev_list}")
        if self.scratch_dev is not None:
            if not is_writable(self.scratch_dev):
                logger.warn(
                    f"Scratch device {self.scratch_dev} is not writable."
                )
                self.scratch_dev = None
            else:
                logger.debug(f"Scratch device is {self.scratch_dev}")
        if not is_writable(self.fallback_dev):
            logger.error(
                f"fallback working device {self.fallback_dev} is not writable!"
            )
            sys.exit(1)
        self.allocatable_devices = self.memdev_list.copy()
        if self.scratch_dev is not None:
            self.allocatable_devices.append(self.scratch_dev)
        self.allocations = {k: 0 for k in self.allocatable_devices}
        self.min_free = {k: free_mb(k) for k in self.allocatable_devices}


@attr.s
class DataMailboxes:
    """Pass data to and from on-disk FIFOs."""

    n_boxes = attr.ib()
    mb_dir_path = attr.ib(default=Path("./mailboxes/"))
    file_extension = attr.ib(default=None)

    def write_headers(self, header):
        """Initialize the mailboxes, writing a free-form header."""
        self.mb_dir_path.mkdir(parents=True, exist_ok=True)
        for i in range(self.n_boxes):
            mb_path = self.path_to_mailbox(i)
            with mb_path.open("w") as fh:
                fh.write(header)

    def write_tsv_headers(self, columns, index_name=None):
        """Initialize the mailboxes, writing a tab-separated header."""
        if index_name is None:
            start = "\t"
        else:
            start = f"{index_name}\t"
        colstring = "\t".join(columns)
        self.mb_dir_path.mkdir(parents=True, exist_ok=True)
        for i in range(self.n_boxes):
            mb_path = self.path_to_mailbox(i)
            with mb_path.open("w") as fh:
                fh.write(f"{start}{colstring}\n")

    @contextlib.contextmanager
    def locked_open_for_write(self, box_no):
        """Acquire a lock on a m."""
        mb_path = self.path_to_mailbox(box_no)
        with mb_path.open("a+") as fd:
            fcntl.flock(fd, fcntl.LOCK_EX)
            yield fd
            fcntl.flock(fd, fcntl.LOCK_UN)

    @contextlib.contextmanager
    def open_then_delete(self, box_no, delete=True):
        """Open a mailbox, then delete when done."""
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
        """Remove the mailbox directory."""
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
        headers = [next(fh).rstrip() for fh in self.fh_list]
        if headers != ([header_value] * self.n_merge):
            print("Error in header values.")
            sys.exit(1)
        self.value_vec = ma.masked_array(
            np.zeros(self.n_merge), mask=np.zeros(self.n_merge),
        ).astype(np.uint32)
        self.payloads = np.array(
            [""] * self.n_merge, dtype=f"<U{MAX_LINE_LEN}"
        )
        self._next_vals(np.arange(self.n_merge))

    def _next_vals(self, index_vec):
        """Return the value of the next iterated value."""
        for index in index_vec:
            try:
                splitline = next(self.fh_list[index]).split("\t")
                self.value_vec[index] = int(splitline[0])
                payload = "\t".join(splitline[1:]).strip()
                if len(payload) > MAX_LINE_LEN:
                    print(f"payload line too long in file index{index}")
                    sys.exit(1)
                self.payloads[index] = payload
            except StopIteration:
                self.value_vec.mask[index] = 1

    def _close_all(self):
        """Close all filehandles."""
        for i in range(self.n_merge):
            self.fh_list[i].close()

    def merge(self, merge_obj):
        """Call merge_obj.merge_func for each merge, return merge_obj.results()."""
        while (~self.value_vec.mask).sum() > 1:
            minimum = np.amin(self.value_vec)
            min_vec = self.value_vec == minimum
            where_min = np.where(min_vec)[0]
            count = len(where_min)
            if count > 1:
                merge_obj.merge_func(
                    minimum,
                    count,
                    ma.masked_array(self.payloads, mask=~min_vec),
                )
            self._next_vals(where_min)
        self._close_all()
        return merge_obj.results()
