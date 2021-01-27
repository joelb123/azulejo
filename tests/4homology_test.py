# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
import sys
from pathlib import Path

# third-party imports
import pytest
import sh

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
def test_homology(datadir_mgr, capsys):
    """Test homology clustering, MSA, and tree building."""
    with capsys.disabled():
        with datadir_mgr.in_tmp_dir(
            inpathlist=INGEST_OUTPUTS,
            save_outputs=True,
            outscope="global",
            excludepaths=["logs/"],
            progressbar=True,
        ):
            args = ["-e", SUBCOMMAND, "glycines"]

            print(f"azulejo {' '.join(args)}")
            try:
                azulejo(
                    args,
                    _out=sys.stderr,
                )
            except sh.ErrorReturnCode as errors:
                print(errors)
                pytest.fail("Homology clustering failed")
            for filestring in find_homology_files():
                assert Path(filestring).stat().st_size > 0
