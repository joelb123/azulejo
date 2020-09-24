# -*- coding: utf-8 -*-
"""Tests for installer function."""
# standard library imports
import sys

# third-party imports
import pytest
import sh

# module imports
from . import help_check
from . import print_docstring
from . import working_directory

# global constants
azulejo = sh.Command("azulejo")
SUBCOMMAND = "install"


def test_subcommand_help():
    """Test subcommand help message."""
    help_check(SUBCOMMAND)


@print_docstring()
def test_installer(tmp_path):
    """Test installer check and build functions."""
    try:
        output = azulejo([SUBCOMMAND])
    except sh.ErrorReturnCode as errors:
        print(errors)
        pytest.fail("Installer test failed")
    assert "muscle" in output
    assert "usearch" in output
    if "All dependencies" in output:
        print("All dependencies already installed.")
    else:
        with working_directory(tmp_path):
            try:
                azulejo([SUBCOMMAND, "-y", "all"], _out=sys.stdout)
            except sh.ErrorReturnCode as errors:
                print(errors)
                pytest.fail("Build failed")
            assert "All dependencies" in azulejo(["install"])
