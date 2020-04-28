# -*- coding: utf-8 -*-
"""
Constants and functions in common across modules
"""
# standard library imports
from pathlib import Path

__NAME__ = "azulejo"
STATFILE_SUFFIX = f"-{__NAME__}_stats.tsv"
ANYFILE_SUFFIX = f"-{__NAME__}_ids-any.tsv"
ALLFILE_SUFFIX = f"-{__NAME__}_ids-all.tsv"
CLUSTFILE_SUFFIX = f"-{__NAME__}_clusts.tsv"
SEQ_FILE_TYPE = "fasta"


def get_paths_from_file(f, must_exist=True):
    inpath = Path(f).expanduser().resolve()
    if must_exist and not inpath.exists():
        raise FileNotFoundError(f)
    dirpath = inpath.parent
    return inpath, dirpath
