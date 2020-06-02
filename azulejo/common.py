# -*- coding: utf-8 -*-
"""Constants and functions in common across modules."""
# standard library imports
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


def fasta_records(filepath):
    """Count the number of records in a FASTA file."""
    count = 0
    next_pos = 0
    angle_bracket = bytes(">", "utf-8")
    with filepath.open("r+b") as fh:
        mm = mmap.mmap(fh.fileno(), 0)
        size = mm.size()
        next_pos = mm.find(angle_bracket, next_pos)
        while next_pos != -1 and next_pos < size:
            count += 1
            next_pos = mm.find(angle_bracket, next_pos + 1)
    return count, size


def fasta_headers(filepath):
    """Return FASTA headers as list of bytes."""
    next_pos = 0
    headers = []
    with filepath.open("r+b") as fh:
        mm = mmap.mmap(fh.fileno(), 0)
        size = mm.size()
        non_header_size = size
        next_pos = mm.find(b">", next_pos)
        while next_pos != -1 and next_pos < size:
            eol_pos = mm.find(b"\n", next_pos)
            if eol_pos == -1:
                break
            header = mm[next_pos + 1 : eol_pos]
            non_header_size -= eol_pos - next_pos
            headers.append(header)
            next_pos = mm.find(b">", eol_pos + 1)
    return headers, non_header_size


def protein_file_stats_filename(setname):
    """Return the name of the protein stat file."""
    return f"{setname}-protein_files.tsv"


def protein_properties_filename(filestem):
    """Return the name of the protein properties file."""
    return f"{filestem}-proteins.tsv"


def homo_degree_dist_filename(filestem):
    """Return the name of the homology degree distribution file."""
    return f"{filestem}-degreedist.tsv"
