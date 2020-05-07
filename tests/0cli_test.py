# -*- coding: utf-8 -*-
# standard library imports
import contextlib
import os
from pathlib import Path

# third-party imports
import pytest
import sh


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
    print("testing basic cli function")
    with working_directory(tmp_path):
        try:
            output = sh.azulejo()
        except sh.ErrorReturnCode as e:
            print(e)
            pytest.fail("Basic cli test failed")
        print(output)
        assert "Usage:" in output
        assert "Options:" in output
        assert "Commands:" in output


def test_version(tmp_path):
    print("testing version")
    with working_directory(tmp_path):
        try:
            output = sh.azulejo(["--version"])
        except sh.ErrorReturnCode as e:
            print(e)
            pytest.fail(e)
        print(output)
        assert "version" in output
