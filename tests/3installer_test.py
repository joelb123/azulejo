# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# third-party imports
import pytest

# first-party imports
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
    """Test installer check function."""
    with working_directory(tmp_path):
        try:
            output = azulejo([SUBCOMMAND])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Installer test failed")
        assert "muscle" in output
        assert "usearch" in output


@print_docstring()
def test_build(tmp_path):
    """Test building dependencies."""
    with working_directory(tmp_path):
        try:
            output = azulejo([SUBCOMMAND, "-f", "-y", "all"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Build failed")
        assert "All dependencies" in azulejo(["install"])
