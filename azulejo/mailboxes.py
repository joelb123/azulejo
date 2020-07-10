# -*- coding: utf-8 -*-
"""Send and receive on-disk messages through names pips with file locking."""
# standard library imports
import contextlib
import fcntl
import os
import sys
from pathlib import Path

# third-party imports
import attr
from loguru import logger

# module imports
from . import cli


@attr.s
class DataMailboxes:

    """Pass data to and from on-disk FIFOs."""

    n_boxes = attr.ib()
    mb_dir_path = attr.ib(default=Path("./mailboxes/"))
    delete_on_exit = attr.ib(default=True)

    def write_headers(self, header):
        """Initialize the mailboxes, optionally writing header."""
        self.mb_dir_path.mkdir(parents=True, exist_ok=True)
        for i in range(self.n_boxes):
            mb_path = self.mb_dir_path / f"{i}"
            with mb_path.open("w") as fh:
                fh.write(header)

    @contextlib.contextmanager
    def locked_open_for_write(self, box_no):
        """Open a mailbox with an advisory file lock."""
        mb_path = self.mb_dir_path / f"{box_no}"
        with mb_path.open("a+") as fd:
            fcntl.flock(fd, fcntl.LOCK_EX)
            yield fd
            fcntl.flock(fd, fcntl.LOCK_UN)

    @contextlib.contextmanager
    def open_then_delete(self, box_no, delete=True):
        """Open a mailbox with an advisory file lock."""
        mb_path = self.mb_dir_path / f"{box_no}"
        with mb_path.open("r") as fd:
            yield fd
            if delete:
                mb_path.unlink()

    def __del__(self, force=False):
        """Remove the mailbox directory. """
        if force or self.delete_on_exit:
            file_list = list(self.mb_dir_path.glob("*"))
            for file in file_list:
                file.unlink()
            self.mb_dir_path.rmdir()
