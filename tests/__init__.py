# -*- coding: utf-8 -*-
"""Base for pytest testing."""
# standard library imports
import contextlib
import os
from pathlib import Path


@contextlib.contextmanager
def working_directory(path):
    """Change working directory in context."""
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
