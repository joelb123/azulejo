# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# third-party imports
import pytest
import sh

from . import working_directory

# global constants
azulejo = sh.Command("azulejo")


def test_installer(tmp_path):
    """Test basic cli function."""
    with working_directory(tmp_path):
        try:
            output = azulejo(["install"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail("Installer test failed")
        print(output)
        assert "muscle" in output
        assert "usearch" in output


def test_build(tmp_path):
    """Test version command."""
    with working_directory(tmp_path):
        try:
            output = azulejo(["install", "all"])
        except sh.ErrorReturnCode as errors:
            print(errors)
            pytest.fail(errors)
        print(output)
        output = azulejo(["install"])
        assert "muscle version at recommended" in output
        assert "usearch version at recommended" in output
