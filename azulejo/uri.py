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
import toml
from requests_download import download as request_download
from loguru import logger
from uri import URI
from xopen import xopen

# global constants
FILE_URI = "file://"
FTP_URI = "ftp://"
REQUIRED_NAMES = ("fasta", "gff")
COMPRESSED_EXTS = ("gz", "bz2", "xz")
FILE_EXTS = {"fasta": ("faa", "fa", "fasta"), "gff": ("gff3", "gff")}
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


def read_input_table(toml_path):
    """Read a TOML of file URIs, return a dicts.

    See the documentation for the characteristics of this
    file.
    """
    try:
        input_dict = toml.load(toml_path)
    except TypeError:
        logger.error(f"Error in filename {toml_path}")
        sys.exit(1)
    except toml.TomlDecodeError:
        logger.error(f"File {toml_path} is not valid TOML.")
        sys.exit(1)
    for key in input_dict:
        for name in REQUIRED_NAMES:
            if name not in input_dict[key]:
                logger.error(
                    f'"{name}" entry not found in file "{toml_path}" entry "{key}"'
                )
                sys.exit(1)
        if "uri" not in input_dict[key]:
            uri = FILE_URI
        else:
            uri = URI(input_dict[key]["uri"])
            del input_dict[key]["uri"]
        transport = str(uri).split(":")[0]
        if transport not in TRANSPORTS:
            logger.error(
                f"Unrecognized transport {transport} in line {idx} URI"
            )
            sys.exit(1)
        for file_type in REQUIRED_NAMES:
            if uri == FILE_URI:
                url = FILE_URI + input_dict[key][file_type]
            else:
                url = str(uri / input_dict[key][file_type])
            input_dict[key][file_type] = url
            if transport == "file":
                file_path = Path(input_dict[key][file_type][len(FILE_URI) :])
                if not file_path.exists():
                    logger.error(
                        f'{file_type} file in {key}, "{file_path}" does not exist'
                    )
                    sys.exit(1)
    return input_dict


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
def read_from_url(url):
    """Read from a URI, transparently decompressing if compressed."""
    filename = url.split("/")[-1]
    if url.startswith(FILE_URI):
        yield xopen(url[len(FILE_URI) :])
    # Not a local file, download to temp directory, cleaning up.
    dirpath = tempfile.mkdtemp()
    logger.debug(f"Downloading {url} via to {dirpath}")

    def cleanup():
        shutil.rmtree(dirpath)

    with _cd(dirpath, cleanup):
        if url.startswith(FTP_URI):
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
                print(f"downloading {url}####################")
                r = request_download(url, filename)
            except:
                logger.error(f"HTTP download error for {url}")
            yield xopen(filename)
