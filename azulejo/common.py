# -*- coding: utf-8 -*-
"""
Constants and functions in common across modules
"""
import sys
from pathlib import Path

__NAME__ = 'azulejo'
STATFILE_SUFFIX = '-%s_stats.tsv' %__NAME__
ANYFILE_SUFFIX = '-%s_ids-any.tsv' %__NAME__
ALLFILE_SUFFIX = '-%s_ids-all.tsv' %__NAME__
CLUSTFILE_SUFFIX = '-%s_clusts.tsv' %__NAME__


def get_paths_from_file(f, must_exist=True):
    inpath = Path(f).expanduser().resolve()
    if must_exist and not inpath.exists():
        raise FileNotFoundError(f)
    dirpath = inpath.parent
    return inpath, dirpath