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
from . import W05_INPUTS
from . import W82_INPUTS
from . import help_check
from . import print_docstring

# global constants
azulejo = sh.Command("azulejo")
TOML_FILE = "glyma+glyso.toml"
SUBCOMMAND = "ingest"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_data_ingestion(datadir_mgr, capsys):
    """Test ingesting sequence data."""
    with capsys.disabled():
        with datadir_mgr.in_tmp_dir(
            inpathlist=W05_INPUTS + W82_INPUTS + [TOML_FILE],
            save_outputs=True,
            outscope="global",
            excludepaths=["logs/"],
        ):
            args = ["-e", SUBCOMMAND, TOML_FILE]
            print(f"azulejo {' '.join(args)}")
            try:
                azulejo(
                    args, _out=sys.stderr,
                )
            except sh.ErrorReturnCode as errors:
                print(errors)
                pytest.fail("Ingestion failed")
            for filestring in INGEST_OUTPUTS:
                assert Path(filestring).exists()
