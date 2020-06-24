# -*- coding: utf-8 -*-
"""Constants and functions in common across modules."""
# standard library imports
import contextlib
import json
import mmap
import os
from pathlib import Path

NAME = "azulejo"
STATFILE_SUFFIX = f"-{NAME}_stats.tsv"
ANYFILE_SUFFIX = f"-{NAME}_ids-any.tsv"
ALLFILE_SUFFIX = f"-{NAME}_ids-all.tsv"
CLUSTFILE_SUFFIX = f"-{NAME}_clusts.tsv"
SEQ_FILE_TYPE = "fasta"

GFF_EXT = "gff3"
FAA_EXT = "faa"
FNA_EXT = "fna"


def cluster_set_name(stem, identity):
    """Get a setname that specifies the %identity value.."""
    if identity == 1.0:
        digits = "10000"
    else:
        digits = (f"{identity:.4f}")[2:]
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
                    f"Header format is bad in {filepath} header {len(properties_dict)+1}"
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
