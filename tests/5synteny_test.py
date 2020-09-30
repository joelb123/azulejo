# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
import time
from pathlib import Path

# third-party imports
import sh

# module imports
from . import HOMOLOGY_OUTPUTS
from . import SYNTENY_OUTPUTS
from . import find_homology_files
from . import help_check
from . import print_docstring
from . import run_azulejo

# global constants
azulejo = sh.Command("azulejo")
SUBCOMMAND = "synteny"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_synteny(datadir_mgr, capsys):
    """Test synteny anchor construction."""
    with capsys.disabled():
        inpathlist = HOMOLOGY_OUTPUTS + find_homology_files(in_tmp_dir=False)
        with datadir_mgr.in_tmp_dir(
            inpathlist=inpathlist,
            progressbar=True,
            save_outputs=True,
            excludepaths=["logs/"],
            outscope="global",
        ):
            run_azulejo(
                ["-e", SUBCOMMAND, "glycines"], "synteny anchor construction"
            )
            print("Checking that output files exist")
            for filestring in SYNTENY_OUTPUTS:
                assert Path(filestring).exists()
