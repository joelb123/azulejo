# -*- coding: utf-8 -*-
"""azulejo -- tile phylogenetic space with subtrees."""
# standard library imports
import locale
import warnings
from pkg_resources import iter_entry_points

# third-party imports
import click
from click_loguru import ClickLoguru
from click_plugins import with_plugins
from loguru import logger

# module imports
from .common import NAME

# global constants
LOG_FILE_RETENTION = 3
__version__ = "0.7.0"

# set locale so grouping works
for localename in ["en_US", "en_US.utf8", "English_United_States"]:
    try:
        locale.setlocale(locale.LC_ALL, localename)
        break
    except locale.Error:
        continue

# set up logging
click_loguru = ClickLoguru(NAME, __version__, retention=LOG_FILE_RETENTION)
# create CLI
@with_plugins(iter_entry_points(NAME + ".cli_plugins"))
@click_loguru.logging_options
@click.group()
@click_loguru.stash_subcommand()
@click.option(
    "-e",
    "--warnings_as_errors",
    is_flag=True,
    show_default=True,
    default=False,
    help="Treat warnings as fatal.",
)
@click.version_option(version=__version__, prog_name=NAME)
def cli(warnings_as_errors, **unused_kwargs):
    """azulejo -- tiling genes in subtrees across phylogenetic space.

    \b
    For more information, see the homepage at https://github.com/legumeinfo/azulejo

    Written by Joel Berendzen <joelb@ncgr.org>.
    Copyright (C) 2020. National Center for Genome Resources. All rights reserved.
    License: BSD-3-Clause
    """
    if warnings_as_errors:
        logger.warning(
            "Runtime warnings (e.g., from pandas) will cause exceptions"
        )
        warnings.filterwarnings("error")


from .analysis import analyze_clusters  # isort:skip
from .analysis import length_std_dist  # isort:skip
from .analysis import outlier_length_dist  # isort:skip
from .core import add_singletons  #  isort:skip
from .core import adjacency_to_graph  #  isort:skip
from .core import cluster_in_steps  #  isort:skip
from .core import clusters_to_histograms  #  isort:skip
from .core import combine_clusters  #  isort:skip
from .core import compare_clusters  #  isort:skip
from .core import prepare_protein_files  #  isort:skip
from .core import usearch_cluster  #  isort:skip
from .synteny import annotate_homology  # isort:skip
from .synteny import join_protein_position_info  # isort:skip
from .synteny import synteny_anchors  # isort:skip
from .synteny import dagchainer_synteny  # isort:skip
from .synteny import proxy_genes  # isort:skip
