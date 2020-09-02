# -*- coding: utf-8 -*-
"""Base for pytest testing."""
# standard library imports
import contextlib
import functools
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


def print_docstring():
    """Decorator to print a docstring."""

    def decorator(func):
        """Define decorator"""

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Print docstring and call function"""
            print(func.__doc__)
            return func(*args, **kwargs)

        return wrapper

    return decorator
