# -*- coding: utf-8 -*-
"""Constants and functions in common across modules."""
# standard library imports
import contextlib
import mmap
import os
import sys
import tempfile
from pathlib import Path

# third-party imports
import numpy as np
import pandas as pd
import xxhash
from loguru import logger as loguru_logger
from memory_tempfile import MemoryTempfile


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

CLUSTER_FILETYPE = "parq"
CLUSTERS_FILE = "clusters.parq"
CLUSTERSYN_FILE = "clusters.syn.parq"
FRAGMENTS_FILE = "fragments.tsv"
ANCHOR_HIST_FILE = "anchor_hist.tsv"
HOMOLOGY_FILE = "proteins.hom.parq"
PROTEOMES_FILE = "proteomes.tsv"
PROTEOMOLOGY_FILE = "proteomes.hom.parq"
PROTEOSYN_FILE = "proteomes.hom.syn.parq"
PROTEINS_FILE = "proteins.parq"
SYNTENY_FILE = "proteins.hom.syn.parq"
ANCHORS_FILE = "anchors.tsv"
SYNTENY_FILETYPE = "tsv"
COLLECTION_FILE = "collection.json"
COLLECTION_HOM_FILE = "collection.hom.json"
COLLECTION_SYN_FILE = "collection.hom.syn.json"
EXTERNAL_CLUSTERS_FILE = "clusters-external.tsv"

# fragment-name defs
PLASTID_STARTS = ["chromop", "chl", "mt", "mi", "rh", "mu", "le", "pl"]
CHROMOSOME_SYNONYMS = ["chromosome", "chrom", "chro", "gs", "gm"]
ALTERNATE_ABBREV = "alt"
CHROMOSOME_ABBREV = "chr"
SCAFFOLD_SYNONYMS = ["scaffold", "scaf", "sca"]
SCAFFOLD_ABBREV = "sc"

# synteny codes
UNAMBIGUOUS_CODE = "U"
DISAMBIGUATED_CODE = "D"
INDIRECT_CODE = "I"
LOCALLY_UNAMBIGUOUS_CODE = "L"
NON_AMBIGUOUS_CODE = "N"
AMBIGUOUS_CODE = "A"
CODE_DICT = {
    UNAMBIGUOUS_CODE: "unambiguous",
    DISAMBIGUATED_CODE: "disambiguated",
    INDIRECT_CODE: "indirectly unambiguous",
    LOCALLY_UNAMBIGUOUS_CODE: "locally unambiguous",
    NON_AMBIGUOUS_CODE: "non-ambiguous",
    AMBIGUOUS_CODE: "ambiguous",
}
DIRECTIONAL_CATEGORY = pd.CategoricalDtype(categories=["-", "+"])
YES_NO = pd.CategoricalDtype(categories=["y", "n"])
SYNTENY_CATEGORY = pd.CategoricalDtype(categories=CODE_DICT.keys())
DEFAULT_DTYPE = pd.UInt32Dtype()
NONDEFAULT_DTYPES = {
    "anchor.subframe.ok": pd.BooleanDtype(),
    "code": SYNTENY_CATEGORY,
    "fasta_url": pd.StringDtype(),
    "gff_url": pd.StringDtype(),
    "frag.direction": DIRECTIONAL_CATEGORY,
    "frag.id": pd.CategoricalDtype(),
    "frag.is_chr": YES_NO,
    "frag.is_plas": YES_NO,
    "frag.is_scaf": YES_NO,
    "frag.len": pd.UInt64Dtype(),
    "frag.orig_id": pd.StringDtype(),
    "frag.start": pd.UInt64Dtype(),
    "gff.feature": pd.CategoricalDtype(),
    "gff.id": pd.CategoricalDtype(),
    "path": pd.CategoricalDtype(),
    "phy.*": pd.CategoricalDtype(),
    "preference": pd.StringDtype(),
    "prot.m_start": pd.BooleanDtype(),
    "prot.no_stop": pd.BooleanDtype(),
    "prot.seq": pd.StringDtype(),
    "syn.anchor.direction": DIRECTIONAL_CATEGORY,
    "syn.code": SYNTENY_CATEGORY,
    "syn.orthogenomic_pct": "float64",
    "syn.hash.footprint": pd.UInt32Dtype(),
    "syn.hash.direction": DIRECTIONAL_CATEGORY,
    "val": "float64",
    # patterns are matched in order after checking for exact matches
    "patterns": [
        {"start": "phy.", "type": pd.CategoricalDtype()},
        {"start": "pct_", "end": "_pct", "type": "float64"},
        {"start": "memb", "type": pd.StringDtype()},
    ],
}

MEGABYTES = 1024.0 * 1024.0

# logger class for use in testing


class PrintLogger:
    """This logger only prints, for testing only."""

    def __init__(self, level):
        try:
            self.level = int(level)
        except ValueError:
            if level.lower() == "debug":
                self.level = 10
            elif level.lower() == "info":
                self.level = 20
            elif level.lower() == "warning":
                self.level = 30
            elif level.lower() == "error":
                self.level = 40
            else:
                self.level = 20

    def debug(self, message):
        if self.level <= 10:
            print(f"Debug: {message}")

    def info(self, message):
        if self.level <= 20:
            print(message)

    def warning(self, message):
        if self.level <= 30:
            print(f"Warning: {message}")

    def error(self, message):
        if self.level <= 40:
            print(f"ERROR: {message}")


def is_writable(dev):
    """Returns whether a device is writable or not."""
    try:
        testdir = tempfile.mkdtemp(prefix=NAME + "-", dir=append_slash(dev))
        Path(testdir).rmdir()
    except OSError:
        return False
    return True


# global variables that are set depending on envvars
# Don't use loguru, just print.  Useful for testing.
if "LOG_TO_PRINT" in os.environ:
    logger = PrintLogger(os.environ["LOG_TO_PRINT"])
else:
    logger = loguru_logger
# Update period on spinners.  Also useful for testing.
if "SPINNER_UPDATE_PERIOD" in os.environ:
    try:
        SPINNER_UPDATE_PERIOD = float(os.environ["SPINNER_UPDATE_PERIOD"])
    except ValueError:
        SPINNER_UPDATE_PERIOD = 5.0
else:
    SPINNER_UPDATE_PERIOD = 1.0
# Fast scratch disk (e.g., SSD or /dev/shm), if other than /tmp
if "SCRATCH_DEV" in os.environ and is_writable(os.environ["SCRATCH_DEV"]):
    SCRATCH_DEV = os.environ["SCRATCH_DEV"]
else:
    SCRATCH_DEV = "/tmp"
# Build disk for installer, needs to allow exe bit set
if "BUILD_DEV" in os.environ and is_writable(os.environ["BUILD_DEV"]):
    BUILD_DEV = os.environ["BUILD_DEV"]
else:
    BUILD_DEV = MemoryTempfile(
        preferred_paths=["/run/user/{uid}"], fallback=True
    ).get_usable_mem_tempdir_paths()[0]


def enforce_canonical_dtypes(frame):
    """Enforce that dtypes of columns meet expectations."""
    for col in frame.columns:
        if col.startswith("tmp."):
            continue
        column_type = frame[col].dtype
        should_be_type = DEFAULT_DTYPE
        if col in NONDEFAULT_DTYPES:
            should_be_type = NONDEFAULT_DTYPES[col]
        else:
            for pattern_dict in NONDEFAULT_DTYPES["patterns"]:
                if "start" in pattern_dict:
                    if col.startswith(pattern_dict["start"]):
                        should_be_type = pattern_dict["type"]
                        break
                if "end" in pattern_dict:
                    if col.endswith(pattern_dict["end"]):
                        should_be_type = pattern_dict["type"]
                        break
        try:
            is_correct_type = column_type == should_be_type
        except TypeError:
            is_correct_type = False
        if not is_correct_type:
            try:
                frame[col] = frame[col].astype(should_be_type)
            except (ValueError, TypeError) as cast_err:
                logger.warning(f"Cannot cast {col} to {should_be_type}")
                logger.warning(cast_err)
    return frame


def free_mb(dev):
    """"Return the number of free MB on dev."""
    fs_stats = os.statvfs(dev)
    free_space_mb = int(
        np.rint((fs_stats.f_bsize * fs_stats.f_bfree) / MEGABYTES)
    )
    return free_space_mb


def append_slash(dev):
    """Append a final slash, if needed."""
    if not dev.endswith("/"):
        dev += "/"
    return dev


def disk_usage_mb(pathstr):
    """Calculate the size used in MB by a path."""
    path = Path(pathstr)
    if not path.exists():
        logger.error(f"Path '{path}' does not exist for disk usage")
        return 0
    total_size = np.array(
        [p.stat().st_size for p in path.rglob("*") if not p.is_symlink()]
    ).sum()
    return int(np.rint(total_size / MEGABYTES))


class MinSpaceTracker:
    """Keep track of the minimum space available on a (memory) device."""

    def __init__(self, device):
        """Remember the device and initialize minimum space."""
        self.device = Path(device)
        self.initial_space = free_mb(self.device)
        self.min_space = self.initial_space

    def check(self):
        """Update the minimimum space available."""
        self.min_space = min(self.min_space, free_mb(self.device))

    def report_min(self):
        """Report the minimum space available."""
        return self.min_space

    def report_used(self):
        """Report change from initial space."""
        return self.initial_space - self.min_space


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
    """Return a dot-separated pathstring as a path."""
    return Path("/".join(dotpath.split(".")))


def fasta_records(filepath):
    """Count the number of records in a FASTA file."""
    count = 0
    next_pos = 0
    angle_bracket = bytes(">", "utf-8")
    memory_map = TrimmableMemoryMap(filepath, access=mmap.ACCESS_READ)
    with memory_map.map() as mem_map:
        size = memory_map.size
        next_pos = mem_map.find(angle_bracket, next_pos)
        while next_pos != -1 and next_pos < size:
            count += 1
            next_pos = mem_map.find(angle_bracket, next_pos + 1)
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


def sort_proteome_frame(frame):
    """Sort a proteome frame by preference and frag.max and renumber."""
    frame = frame.copy()
    if frame.index.name == "path":
        frame["path"] = frame.index
    frame.sort_values(
        by=["preference", "frag.max"], ascending=[True, False], inplace=True
    )
    frame["order"] = range(len(frame))
    frame.set_index("order", inplace=True)
    return frame


def remove_tmp_columns(frame):
    """Remove any columns in a data frame that begin with 'tmp.'."""
    drop_cols = [col for col in frame.columns if col.startswith("tmp.")]
    if len(drop_cols) != 0:
        return frame.drop(drop_cols, axis=1)
    return frame


def write_tsv_or_parquet(
    frame,
    filepath,
    compression=DEFAULT_PARQUET_COMPRESSION,
    float_format="%.2f",
    desc=None,
    remove_tmp=True,
    sort_cols=True,
    enforce_types=True,
):
    """Write either a TSV or a parquet file by file extension."""
    filepath = Path(filepath)
    ext = filepath.suffix.lstrip(".")
    if desc is not None:
        file_desc = f"{desc} file"
        logger.debug(f'Writing {file_desc} "{filepath}')
    if remove_tmp:
        frame = remove_tmp_columns(frame)
    if enforce_types:
        frame = enforce_canonical_dtypes(frame)
    if sort_cols:
        frame = frame[sorted(frame.columns)]
    if ext in PARQUET_EXTENSIONS:
        frame.to_parquet(filepath, compression=compression)
    elif ext in TSV_EXTENSIONS:
        frame.to_csv(filepath, sep="\t", float_format=float_format)
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
    if ext in TSV_EXTENSIONS:
        frame = pd.read_csv(filepath, sep="\t", index_col=0).convert_dtypes()
        return enforce_canonical_dtypes(frame)
    logger.error(f"Unrecognized file extensions {ext} in {filepath}")
    sys.exit(1)


def log_and_add_to_stats(stats, new_stats):
    """Print stats info and write to stats file."""
    with pd.option_context(
        "display.max_rows",
        None,
        "display.max_columns",
        None,
        "display.float_format",
        "{:,.1f}%".format,
    ):
        logger.info(new_stats)
    overlap_cols = list(set(stats.columns) & set(new_stats.columns))
    return pd.concat([stats.drop(columns=overlap_cols), new_stats], axis=1)


def bool_to_y_or_n(bool_arr):
    """Convert boolean to T/F value"""
    boolean_dict = {True: "y", False: "n"}
    bool_ser = pd.Series(bool_arr)
    return bool_ser.map(boolean_dict)


def y_or_n_to_bool(bool_ser):
    """Convert boolean to T/F value"""
    tf_dict = {"y": True, False: "n"}
    return bool_ser.map(tf_dict).astype(bool)


def hash_array(kmer):
    """Return a hash of a numpy array."""
    return xxhash.xxh32_intdigest(kmer.tobytes())


def calculate_adjacency_group(index_series, frag_series):
    """Calculate an adjacency group numger."""
    index_fr = pd.DataFrame({"index": index_series, "fragment": frag_series})
    n_prot = len(index_fr)
    adj_gr_count = 0
    was_adj = False
    index_fr["i"] = range(n_prot)
    adj_group = np.array([np.nan] * n_prot)
    for unused_group, subframe in index_fr.groupby(by=["fragment"]):
        if len(subframe) == 1:
            continue
        last_pos = -2
        last_row = None
        if was_adj:
            adj_gr_count += 1
        was_adj = False
        for unused_i, row in subframe.iterrows():
            row_no = row["i"]
            if row["index"] == last_pos + 1:
                if not was_adj:
                    adj_group[last_row] = adj_gr_count
                was_adj = True
                adj_group[row_no] = adj_gr_count
            else:
                if was_adj:
                    adj_gr_count += 1
                    was_adj = False
            last_pos = row["index"]
            last_row = row_no
    if was_adj:
        adj_gr_count += 1
    adj_arr = pd.Series(
        adj_group, dtype=pd.UInt32Dtype(), index=index_series.index
    )
    n_adj = n_prot - adj_arr.isnull().sum()
    return n_adj, adj_gr_count, adj_arr
