# -*- coding: utf-8 -*-
"""Constants and functions in common across modules."""
# standard library imports
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


def get_paths_from_file(filepath, must_exist=True):
    """Given a string filepath,, return the resolved path and parent."""
    inpath = Path(filepath).expanduser().resolve()
    if must_exist and not inpath.exists():
        raise FileNotFoundError(filepath)
    dirpath = inpath.parent
    return inpath, dirpath


def protein_file_stats_filename(setname):
    """Return a string that is the name of the protein stat file."""
    return f"{setname}-protein_file_stats.tsv"


def protein_properties_filename(filestem):
    """Return a string that is the name of the protein properties file."""
    return f"{filestem}-protein_stats.tsv"
