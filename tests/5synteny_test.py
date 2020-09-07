# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
import sys
from pathlib import Path

# third-party imports
import pytest
import sh

# module imports
from . import HOMOLOGY_OUTPUTS
from . import SYNTENY_OUTPUTS
from . import find_homology_files
from . import help_check
from . import print_docstring

# global constants
azulejo = sh.Command("azulejo")
SUBCOMMAND = "synteny-anchors"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_synteny_anchors(datadir_mgr):
    """Test synteny anchor construction."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=HOMOLOGY_OUTPUTS + find_homology_files(in_tmp_dir=False),
        save_outputs=True,
        outscope="global",
        excludepaths=["logs/"],
    ):
        try:
            output = azulejo(["-e", "-q", SUBCOMMAND, "glycines"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Synteny anchor construction failed")
        print(f"output={output}")
        for filestring in SYNTENY_OUTPUTS:
            assert Path(filestring).exists()
