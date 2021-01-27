# -*- coding: utf-8 -*-
"""Base for pytest testing."""
# standard library imports
import contextlib
import functools
import os
import sys
from pathlib import Path

# third-party imports
import pytest
import sh
from sh import ErrorReturnCode

# global constants
W05_INPUTS = [
    "glyso.W05.gnm1.ann1.T47J.protein_primaryTranscript.faa",
    "glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3",
]

W82_INPUTS = [
    "glyma.Wm82.gnm4.ann1.T8TQ.gene_models_main.gff3",
    "glyma.Wm82.gnm4.ann1.T8TQ.protein_primaryTranscript.faa",
]

SET_DIR = "glycines"
PROT_SUBDIRS = ["glyso/W05/", "glyma/"]
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
            "homology_clusters.parq",
            "proteomes.hom.parq",
            "homology_cluster_hist.tsv",
            "homology-stats.tsv",
        ]
        + [f"{subdir}proteins.hom.parq" for subdir in PROT_SUBDIRS]
    )
]

SYNTENY_OUTPUTS = [
    f"{SET_DIR}/{f}"
    for f in (
        ["proteomes.hom.syn.parq"]  # "clusters.syn.parq",
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
        output = sh.azulejo(help_command)
    except ErrorReturnCode as errors:
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


def run_azulejo(args, component):
    """Run azulejo with args."""
    command_string = " ".join(args)
    print(f"Testing {component} with" + f'"azulejo {command_string}"')
    try:
        sh.azulejo(
            args,
            _out=sys.stderr,
        )
    except ErrorReturnCode as errors:
        print(errors)
        pytest.fail(f"{component} failed")
