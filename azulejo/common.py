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
DEFAULT_PARQUET_COMPRESSION = "ZSTD"
# Sizes and minimum read times with various compressions
# for a file with one proteome on a system with M.2 SSD disk
# under pyarrow 1.0.0 into pandas 1.1.0:
#   "NONE": 43MB, 1.8s
#   "ZSTD": 13M, 1.8s
# "SNAPPY": 29 MB, 1.8s
# "BROTLI": 13 MB, 1.9s
#    "LZ4": 23MB, (disabled under pyarrow 1.0.0, was about like brotli under 0.17)
#   "GZIP": 14 MB, 2.1 s
#    "LZO": not supported
#    "BZ2": not supported
# In addition, the ingest process took 28.8s with None, and
# 28.4 s with ZSTD, probably due to writing less data.
# With its 70% compression factor,  ZSTD can be expected to
# perform even better relative to uncompressed and snappy
# on production systems with slower disks for which
# cache is not warmed up (as mine was in this test).
# So ZSTD seems a clear choice for now.

PARQUET_EXTENSIONS = ["parquet", "pq", "parq"]
TSV_EXTENSIONS = ["tsv"]
SAVED_INPUT_FILE = "input.toml"

# Changing the extension of these files will change the type of file written.
# TSV files, though readable/editable, do not give the written values back.
# Parquet is also ~100X faster.
4
CLUSTER_FILETYPE = "parq"
CLUSTERS_FILE = "clusters.parq"
FRAGMENTS_FILE = "fragments.tsv"
ANCHOR_HIST_FILE = "anchor_hist.tsv"
HOMOLOGY_FILE = "proteins.hom.parq"
PROTEOMES_FILE = "proteomes.tsv"
PROTEOMOLOGY_FILE = "proteomes.hom.parq"
PROTEOSYN_FILE = "proteomes.hom.syn.parq"
PROTEINS_FILE = "proteins.parq"
UNRENAMED_PROTEINS_FILE = "proteins-unrenamed.parq"
SYNTENY_FILE = "proteins.hom.syn.parq"

# fragment-name defs
PLASTID_STARTS = ["chromop", "chl", "mt", "mi", "rh", "mu", "le", "pl"]
CHROMOSOME_SYNONYMS = ["chromosome", "chrom", "chro", "gs", "gm"]
ALTERNATE_ABBREV = "alt"
CHROMOSOME_ABBREV = "chr"
SCAFFOLD_SYNONYMS = ["scaffold", "scaf", "sca"]
SCAFFOLD_ABBREV = "sc"

DIRECTIONAL_CATEGORY = pd.CategoricalDtype(categories=["-", "+"])
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
                self.filehandle.close()
            else:
                self.map_obj.close()
                self.filehandle.close()


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
    df = df.copy()
    if df.index.name == "path":
        df["path"] = df.index
    df.sort_values(
        by=["preference", "frag.max"], ascending=[True, False], inplace=True
    )
    df["order"] = range(len(df))
    df.set_index("order", inplace=True)
    return df


def remove_tmp_columns(df):
    """Remove any columns in a data frame that begin with 'tmp.'."""
    drop_cols = [col for col in df.columns if col.startswith("tmp.")]
    if len(drop_cols) != 0:
        return df.drop(drop_cols, axis=1)
    return df


def write_tsv_or_parquet(
    df,
    filepath,
    compression=DEFAULT_PARQUET_COMPRESSION,
    float_format="%.2f",
    desc=None,
    remove_tmp=True,
):
    """Write either a TSV or a parquet file by file extension."""
    filepath = Path(filepath)
    ext = filepath.suffix.lstrip(".")
    if desc is not None:
        file_desc = f"{desc} file"
        logger.debug(f'Writing {file_desc} "{filepath}')
    if remove_tmp:
        df = remove_tmp_columns(df)
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
        return pd.read_parquet(filepath).convert_dtypes()
    elif ext in TSV_EXTENSIONS:
        df = pd.read_csv(filepath, sep="\t", index_col=0)
        return df.convert_dtypes()
    else:
        logger.error(f"Unrecognized file extensions {ext} in {filepath}")
        sys.exit(1)
