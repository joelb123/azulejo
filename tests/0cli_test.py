# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# third-party imports
import pytest
import sh

# module imports
from . import print_docstring
from . import working_directory

# global constants
azulejo = sh.Command("azulejo")


@print_docstring()
def test_cli(tmp_path):
    """Test basic cli function."""
    with working_directory(tmp_path):
        try:
            output = azulejo()
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Basic cli test failed")
        assert "Usage:" in output
        assert "Options:" in output
        assert "Commands:" in output


@print_docstring()
def test_version(tmp_path):
    """Test version command."""
    with working_directory(tmp_path):
        try:
            output = azulejo(["--version"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(errors)
        assert "version" in output


@print_docstring()
def test_taxonomy(tmp_path):
    """Test taxonomy rank check command."""
    with working_directory(tmp_path):
        try:
            output = azulejo(["check-taxonomic-rank"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(errors)
        assert int(azulejo(["check-taxonomic-rank", "species"])) == 130
        assert int(azulejo(["check-taxonomic-rank", "subspecies"])) == 131
        assert int(azulejo(["check-taxonomic-rank", "superspecies"])) == 128
