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
from . import INGEST_OUTPUTS
from . import help_check
from . import print_docstring

# global constants
azulejo = sh.Command("azulejo")
SUBCOMMAND = "cluster-build-trees"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_cluster_build_trees(datadir_mgr):
    """Test homology clustering, MSA, and tree building."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=INGEST_OUTPUTS,
        save_outputs=True,
        outscope="global",
        excludepaths=["logs/"],
    ):
        try:
            output = azulejo(["-e", SUBCOMMAND, "glycines"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Homology clustering failed")
        print(f"output={output}")
        for filestring in HOMOLOGY_OUTPUTS:
            if not Path(filestring).exists():
                print(f"{filestring} not found")
