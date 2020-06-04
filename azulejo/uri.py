# -*- coding: utf-8 -*-
"""Constants and functions in common across modules."""
# standard library imports\
import contextlib
import os
import sys
import shutil
import tempfile
import urllib.request as request
from urllib.error import URLError
from pathlib import Path

# third-party imports
import pandas as pd
import xopen
from requests_download import download as request_download
from loguru import logger
from uri import URI

# global constants
FILE_URI = "file://"
FTP_URI = "ftp://"
COLUMN_NAMES = ("label", "feature", "protein", "URI")
COMPRESSED_EXTS = ("gz", "bz2", "xz")
FILE_EXTS = {"protein": ("faa", "fa", "fasta"), "feature": ("gff3", "gff")}
TRANSPORTS = ("file", "http", "https", "ftp")


def _compression_type(namedict, filetype):
    """Return True if final extension is compressed type."""
    dot_splits = namedict[filetype].split(".")
    if dot_splits[-1] in COMPRESSED_EXTS:
        _check_ext(dot_splits[-2], filetype)
        return dot_splits[-1]
    _check_ext(dot_splits[-1], filetype)
    return None


def _check_ext(ext, filetype):
    """Check if extension is valid for type."""
    if ext not in FILE_EXTS[filetype]:
        logger.error(f"{ext} is not a valid file extension")
        sys.exit(1)


def read_input_table(table_path):
    """Read a TSV table of file URIs, return a list of dicts.

    See the documentation for the characteristics of this
    file.
    """
    input_table = pd.read_csv(table_path, sep="\t")
    for column_name in COLUMN_NAMES:
        if column_name not in input_table.columns:
            logger.error(
                f'Column name "{column_name}" not found in file "{table_path}'
            )
            sys.exit(1)
    output_list = []
    for idx, row in input_table.iterrows():
        if row["URI"] == "":
            uri = URI(FILE_URI)
        else:
            uri = URI(row["URI"])
        transport = str(uri).split(":")[0]
        if transport not in TRANSPORTS:
            logger.error(
                f"Unrecognized transport {transport} in line {idx} URI"
            )
            sys.exit(1)
        entry_dict = {"label": row["label"]}
        for file_type in ["feature", "protein"]:
            entry_dict[file_type] = str(uri / row[file_type])
            entry_dict[file_type + "_compression"] = _compression_type(
                row, filetype
            )
            if transport == "file":
                file_path = Path(entry_dict[file_type][8:])
                if not file_path.exists():
                    logger.error(
                        f'{file_type} file in row {idx}, "{file_path}" does not exist'
                    )
                    sys.exit(1)
        output_list.append(entry_dict)
    return output_list


@contextlib.contextmanager
def _cd(newdir, cleanup=lambda: True):
    "Change directory with cleanup."
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()


@contextlib.contextmanager
def read_from_uri(url):
    """Read from a URI, transparently decompressing if compressed."""
    filename = url.split("/")[-1]
    if url.beginswith(FILE_URI):
        yield xopen(url[len(FILE_URI) :])
    # Not a local file, download to temp directory, cleaning up.
    dirpath = tempfile.mkdtemp()
    logger.debug(f"Downloading {url} via to {dirpath}")

    def cleanup():
        shutil.rmtree(dirpath)

    with _cd(dirpath, cleanup):
        if url.beginswith(FTP_URI):
            try:
                with contextlib.closing(request.urlopen(url)) as r:
                    with open(filename, "wb") as f:
                        shutil.copyfileobj(r, f)
            except URLError as e:
                if e.reason.find("No such file or directory") >= 0:
                    logger.error(f"FileNotFound from {url}")
                    sys.exit(1)
                else:
                    logger.error(
                        f"FTP download of {url} failed with reason {e.reason}"
                    )
                    sys.exit(1)
                yield xopen(filename)
        else:
            try:
                request_download(url, filename)
            except:
                logger.error(f"HTTP download error for {url}")
            yield xopen(filename)
