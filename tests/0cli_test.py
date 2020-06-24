# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# third-party imports
import pytest
import sh

from . import working_directory

# global constants
azulejo = sh.Command("azulejo")


def test_cli(tmp_path):
    """Test basic cli function."""
    with working_directory(tmp_path):
        try:
            output = azulejo()
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Basic cli test failed")
        print(output)
        assert "Usage:" in output
        assert "Options:" in output
        assert "Commands:" in output


def test_version(tmp_path):
    """Test version command."""
    with working_directory(tmp_path):
        try:
            output = azulejo(["--version"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(errors)
        print(output)
        assert "version" in output


def test_taxonomy(tmp_path):
    """Test version command."""
    with working_directory(tmp_path):
        try:
            output = azulejo(["check-taxonomic-rank"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(errors)
        print(output)
        assert int(azulejo(["check-taxonomic-rank", "species"])) == 130
        assert int(azulejo(["check-taxonomic-rank", "subspecies"])) == 131
        assert int(azulejo(["check-taxonomic-rank", "superspecies"])) == 128
