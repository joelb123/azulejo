# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# standard library imports
import contextlib
import os
from pathlib import Path

# third-party imports
import pytest
import sh

# global constant
azulejo = sh.Command("azulejo")


@contextlib.contextmanager
def working_directory(path):
    """Change working directory in context."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)


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
