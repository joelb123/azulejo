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
def test_install_all(tmp_path, capsys):
    """Test install all required binaries."""
    with capsys.disabled():
        with working_directory(tmp_path):
            try:
                azulejo([SUBCOMMAND, "-y", "all"], _out=sys.stdout)
            except sh.ErrorReturnCode as errors:
                print(errors)
                pytest.fail("Build failed")
            assert "All required" in azulejo(["install"])


@print_docstring()
def test_install_dagchainer_tool(tmp_path, capsys):
    """Test install dagchainer-tool."""
    with capsys.disabled():
        with working_directory(tmp_path):
            try:
                azulejo([SUBCOMMAND, "dagchainer-tool"], _out=sys.stdout)
            except sh.ErrorReturnCode as errors:
                print(errors)
                pytest.fail("Build failed")
            assert "dagchainer-tool.sh at recommended" in azulejo(["install"])


#@print_docstring()
#def test_install_blast(tmp_path, capsys):
#    """Test install BLAST with long ids."""
#    with capsys.disabled():
#        with working_directory(tmp_path):
#            try:
#                azulejo([SUBCOMMAND, "blast-longids"], _out=sys.stdout)
#            except sh.ErrorReturnCode as errors:
#                print(errors)
#                pytest.fail("Build failed")
#            assert "blastn at recommended" in azulejo(["install"])
