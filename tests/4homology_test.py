# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
import sh
from pathlib import Path

# third-party imports
import pytest
import sys

# module imports
from . import INGEST_OUTPUTS
from . import find_homology_files
from . import help_check
from . import print_docstring

# global constants
azulejo = sh.Command("azulejo")
SUBCOMMAND = "homology"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_homology(datadir_mgr, caplog):
    """Test homology clustering, MSA, and tree building."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=INGEST_OUTPUTS,
        save_outputs=True,
        outscope="global",
        excludepaths=["logs/"],
    ):
        try:
            output = azulejo(
                ["-e", SUBCOMMAND, "glycines"], _out=sys.stdout, _tty_out=True
            )
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Homology clustering failed")
        print(output)
        print(caplog.text)
        for filestring in find_homology_files():
            assert Path(filestring).stat().st_size > 0
