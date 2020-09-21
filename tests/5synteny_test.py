# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
import sys
from pathlib import Path

# third-party imports
import pytest
import sh
from click.testing import CliRunner

# first-party imports
from azulejo import cli

# module imports
from . import HOMOLOGY_OUTPUTS
from . import SYNTENY_OUTPUTS
from . import find_homology_files
from . import help_check
from . import print_docstring

# global constants
azulejo = sh.Command("azulejo")
SUBCOMMAND = "synteny"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_synteny(datadir_mgr, caplog):
    """Test synteny anchor construction."""
    runner = CliRunner()
    print(f"caplog={caplog}")
    with datadir_mgr.in_tmp_dir(
        inpathlist=HOMOLOGY_OUTPUTS + find_homology_files(in_tmp_dir=False),
        save_outputs=True,
        outscope="global",
        excludepaths=["logs/"],
    ):
        print("running test")
        try:
            azulejo(
                ["-e", SUBCOMMAND, "glycines"], _out=sys.stdout, _tty_out=True
            )
            # result = runner.invoke(cli, ["-e", SUBCOMMAND, "glycines"])

        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Synteny anchor construction failed")
        # assert result.exit_code == 0
        # print(result.output)
        print(f"checking that output files exist")
        for filestring in SYNTENY_OUTPUTS:
            assert Path(filestring).exists()
        print("done with test, copying output files")
