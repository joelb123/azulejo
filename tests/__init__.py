# -*- coding: utf-8 -*-
"""Base for pytest testing."""
# standard library imports
import contextlib
import functools
import os
from pathlib import Path

# third-party imports
import pytest
import sh

# global constants
W05_INPUTS = [
    "glyso.W05.gnm1.ann1.T47J.protein_primaryTranscript.faa",
    "glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3",
]

SET_DIR = "glycines"
PROT_SUBDIRS = ["glyso/W05/", "glyso/PI483463/", "glyma/"]
N_CLUSTERS = 14617
N_ANCHORS = 32335
INGEST_OUTPUTS = [
    f"{SET_DIR}/{f}"
    for f in (
        ["fragments.tsv", "proteomes.tsv"]
        + [f"{subdir}proteins.parq" for subdir in PROT_SUBDIRS]
    )
]

HOMOLOGY_OUTPUTS = [
    f"{SET_DIR}/{f}"
    for f in (
        [
            "clusters.parq",
            "proteomes.hom.parq",
            "groupkeys-2.tsv",
            "cluster_hist.tsv",
            "homology-stats.tsv",
        ]
        + [f"{subdir}proteins.hom.parq" for subdir in PROT_SUBDIRS]
    )
]

SYNTENY_OUTPUTS = [  # [f"{SET_DIR}/synteny/{i}.parq" for i in range(N_ANCHORS)] +
    f"{SET_DIR}/{f}"
    for f in (
        ["clusters.syn.parq", "proteomes.hom.syn.parq"]
        + [f"{subdir}proteins.hom.syn.parq" for subdir in PROT_SUBDIRS]
    )
]

TSV_TEST_FILE = f"{SET_DIR}/{PROT_SUBDIRS[0]}/proteins.hom.syn.parq"
TSV_OUTPUT_FILE = "proteins.hom.syn.tsv"


def find_homology_files(in_tmp_dir=True):
    """Find the set of homology files, which can vary."""
    if in_tmp_dir:
        prefix = Path(".")
    else:
        prefix = Path("tests/data")
    subdir = f"{SET_DIR}/homology/"
    homology_files = (prefix / subdir).glob("*.parq")
    return [subdir + f.name for f in homology_files]


@contextlib.contextmanager
def working_directory(path):
    """Change working directory in context."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


def help_check(subcommand):
    """Test help function for subcommand."""
    print(f"Test {subcommand} help.")
    if subcommand == "global":
        help_command = ["--help"]
    else:
        help_command = [subcommand, "--help"]
    try:
        output = sh.Command("azulejo")(help_command)
    except sh.ErrorReturnCode as errors:
        print(errors)
        pytest.fail(f"{subcommand} help test failed")
    print(output)
    assert "Usage:" in output
    assert "Options:" in output


def print_docstring():
    """Decorator to print a docstring."""

    def decorator(func):
        """Define decorator"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Print docstring and call function"""
            print(func.__doc__)
            return func(*args, **kwargs)

        return wrapper

    return decorator
