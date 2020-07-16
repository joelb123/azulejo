# -*- coding: utf-8 -*-
"""Constants and functions in common across modules."""
# standard library imports
import contextlib
import json
import mmap
import sys
from pathlib import Path

# third-party imports
import pandas as pd
from loguru import logger

# global constants
NAME = "azulejo"
DEFAULT_PARQUET_COMPRESSION = None
PARQUET_EXTENSIONS = ["parquet", "pq", "parq"]
TSV_EXTENSIONS = ["tsv"]
SAVED_INPUT_FILE = "input.toml"

# Changing the extension of these files will change the type of file written.
# TSV files, though readable/editable, do not give the written values back.
# Parquet is also ~100X faster.
4
CLUSTER_FILETYPE = "tsv"
CLUSTERS_FILE = "clusters.parq"
FRAGMENTS_FILE = "fragments.tsv"
HASH_HIST_FILE = "hash_hist.tsv"
HOMOLOGY_FILE = "proteins+homology.parq"
MERGED_HASH_FILE = "merged_hashes.parq"
PROTEOMES_FILE = "proteomes.tsv"
PROTEOMOLOGY_FILE = "proteomes+homology.parq"
PROTEINS_FILE = "proteins.parq"
SYNTENY_FILE = "proteins+homology+synteny.parq"

# shared functions


def cluster_set_name(stem, identity):
    """Get a setname that specifies the %identity value.."""
    if identity == 1.0:
        digits = "10000"
    else:
        digits = f"{identity:.4f}"[2:]
    return f"{stem}-nr-{digits}"


def get_paths_from_file(filepath, must_exist=True):
    """Given a string filepath,, return the resolved path and parent."""
    inpath = Path(filepath).expanduser().resolve()
    if must_exist and not inpath.exists():
        raise FileNotFoundError(filepath)
    dirpath = inpath.parent
    return inpath, dirpath


class TrimmableMemoryMap:
    """A memory-mapped file that can be resized at the end."""

    def __init__(self, filepath, access=mmap.ACCESS_WRITE):
        """Open the memory-mapped file."""
        self.orig_size = None
        self.size = None
        self.map_obj = None
        self.access = access
        self.filehandle = open(filepath, "r+b")

    def trim(self, start, end):
        """Trim the memory map and mark the nex size."""
        self.map_obj.move(start, end, self.orig_size - end)
        self.size -= end - start
        return self.size

    @contextlib.contextmanager
    def map(self):
        """Open a memory-mapped view of filepath."""
        try:
            self.map_obj = mmap.mmap(
                self.filehandle.fileno(), 0, access=self.access
            )
            self.orig_size = self.map_obj.size()
            self.size = self.orig_size
            yield self.map_obj
        finally:
            if self.access == mmap.ACCESS_WRITE:
                self.map_obj.flush()
                self.map_obj.close()
                self.filehandle.truncate(self.size)


def dotpath_to_path(dotpath):
    "Return a dot-separated pathstring as a path."
    return Path("/".join(dotpath.split(".")))


def fasta_records(filepath):
    """Count the number of records in a FASTA file."""
    count = 0
    next_pos = 0
    angle_bracket = bytes(">", "utf-8")
    memory_map = TrimmableMemoryMap(filepath, access=mmap.ACCESS_READ)
    with memory_map.map() as mm:
        size = memory_map.size
        next_pos = mm.find(angle_bracket, next_pos)
        while next_pos != -1 and next_pos < size:
            count += 1
            next_pos = mm.find(angle_bracket, next_pos + 1)
    return count, size


def parse_cluster_fasta(filepath, trim_dict=True):
    """Return FASTA headers as a dictionary of properties."""
    next_pos = 0
    properties_dict = {}
    memory_map = TrimmableMemoryMap(filepath)
    with memory_map.map() as mm:
        size = memory_map.size
        next_pos = mm.find(b">", next_pos)
        while next_pos != -1 and next_pos < size:
            eol_pos = mm.find(b"\n", next_pos)
            if eol_pos == -1:
                break
            space_pos = mm.find(b" ", next_pos + 1, eol_pos)
            if space_pos == -1:
                raise ValueError(
                    f"Header format is bad in {filepath} header"
                    f" {len(properties_dict)+1}"
                )
            id = mm[next_pos + 1 : space_pos].decode("utf-8")
            payload = json.loads(mm[space_pos + 1 : eol_pos])
            properties_dict[id] = payload
            if trim_dict:
                size = memory_map.trim(space_pos, eol_pos)
            next_pos = mm.find(b">", space_pos)
    return properties_dict


def protein_file_stats_filename(setname):
    """Return the name of the protein stat file."""
    if setname is None:
        return "protein_files.tsv"
    return f"{setname}-protein_files.tsv"


def protein_properties_filename(filestem):
    """Return the name of the protein properties file."""
    if filestem is None:
        return "proteins.tsv"
    return f"{filestem}-proteins.tsv"


def homo_degree_dist_filename(filestem):
    """Return the name of the homology degree distribution file."""
    return f"{filestem}-degreedist.tsv"


def group_key_filename(members):
    """Return the name of the group key file."""
    return f"groupkeys-{members}.tsv"


def sort_proteome_frame(df):
    """Sort a proteome frame by preference and frag.max and renumber."""
    if df.index.name == "path":
        df["path"] = df.index
    df.sort_values(
        by=["preference", "frag.max"], ascending=[True, False], inplace=True
    )
    df["order"] = range(len(df))
    df.set_index("order", inplace=True)
    return df


def write_tsv_or_parquet(
    df,
    filepath,
    compression=DEFAULT_PARQUET_COMPRESSION,
    float_format=None,
    desc=None,
):
    """Write either a TSV or a parquet file by file extension."""
    filepath = Path(filepath)
    ext = filepath.suffix.lstrip(".")
    if desc is not None:
        file_desc = f"{desc} file"
        logger.debug(f'Writing {file_desc} "{filepath}')
    if ext in PARQUET_EXTENSIONS:
        df.to_parquet(filepath, compression=compression)
    elif ext in TSV_EXTENSIONS:
        df.to_csv(filepath, sep="\t", float_format=float_format)
    else:
        logger.error(f"Unrecognized file extension {ext} in {filepath}")
        sys.exit(1)


def read_tsv_or_parquet(filepath):
    """Read either a TSV or a parquet file by file extension."""
    filepath = Path(filepath)
    if not filepath.exists():
        logger.error(f'File "{filepath}" does not exist.')
        sys.exit(1)
    ext = filepath.suffix.lstrip(".")
    if ext in PARQUET_EXTENSIONS:
        return pd.read_parquet(filepath)
    elif ext in TSV_EXTENSIONS:
        return pd.read_csv(filepath, sep="\t", index_col=0)
    else:
        logger.error(f"Unrecognized file extensions {ext} in {filepath}")
        sys.exit(1)
