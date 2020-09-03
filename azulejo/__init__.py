# -*- coding: utf-8 -*-
"""azulejo -- tile phylogenetic space with subtrees."""
# standard library imports
import locale
import os
import warnings
from pathlib import Path
from pkg_resources import iter_entry_points

# third-party imports
import click
from click_plugins import with_plugins
from click_loguru import ClickLoguru

# module imports
from .common import NAME
from .installer import DependencyInstaller

# global constants
LOG_FILE_RETENTION = 3
__version__ = "0.9.8"
INSTALL_ENVIRON_VAR = (  # installs go into "/bin" and other subdirs of this directory
    NAME.upper() + "_INSTALL_DIR"
)
if INSTALL_ENVIRON_VAR in os.environ:
    INSTALL_PATH = Path(os.environ[INSTALL_ENVIRON_VAR])
else:
    INSTALL_PATH = None

MUSCLE_VER = "3.8.1551"
USEARCH_VER = "11.0.667"

# helper functions
def _muscle_version_parser(ver_str):
    """Parse version out of muscle version string."""
    return ver_str.split()[1]


def _usearch_version_parser(ver_str):
    """Parse version out of usearch version string."""
    return ver_str.split()[1].split("_")[0]


DEPENDENCY_DICT = {
    "muscle": {
        "binaries": ["muscle"],
        "tarball": (
            "https://www.drive5.com/muscle/muscle_src_"
            + MUSCLE_VER
            + ".tar.gz"
        ),
        "dir": ".",
        "version": MUSCLE_VER,
        "version_command": ["-version"],
        "version_parser": _muscle_version_parser,
        "make": ["muscle-pgo",],
        "copy_binaries": ["muscle"],
        "license": "public domain",
        "license_restrictive": False,
    },
    "usearch": {
        "binaries": ["usearch"],
        "dir": ".",
        "version": USEARCH_VER,
        "version_command": ["-version"],
        "version_parser": _usearch_version_parser,
        "download_binaries": {
            "linux": "https://www.drive5.com/downloads/usearch11.0.667_i86linux32.gz",
            "macos": "https://www.drive5.com/downloads/usearch11.0.667_i86osx32.gz",
        },
        "make": ["usearch"],
        "copy_binaries": ["usearch"],
        "license": "proprietary non-commercial",
        "license_restrictive": True,
        "license_file": "LICENSE.txt",
    },
}

# set locale so grouping works
for localename in ["en_US", "en_US.utf8", "English_United_States"]:
    try:
        locale.setlocale(locale.LC_ALL, localename)
        break
    except locale.Error:
        continue

# set up logging
click_loguru = ClickLoguru(
    NAME, __version__, retention=LOG_FILE_RETENTION, timer_log_level="info"
)
# create CLI
@with_plugins(iter_entry_points(NAME + ".cli_plugins"))
@click_loguru.logging_options
@click.group(
    epilog="Dependencies: "
    + DependencyInstaller(
        DEPENDENCY_DICT, pkg_name=NAME, install_path=INSTALL_PATH
    ).status()
)
@click_loguru.stash_subcommand()
@click.option(
    "-e",
    "--warnings_as_errors",
    is_flag=True,
    show_default=True,
    default=False,
    help="Treat warnings as fatal.",
    callback=click_loguru.user_global_options_callback,
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    default=True,
    show_default=True,
    help="Process in parallel.",
    callback=click_loguru.user_global_options_callback,
)
@click.version_option(version=__version__, prog_name=NAME)
def cli(warnings_as_errors, parallel, **unused_kwargs):
    """Azulejo -- tiling genes in subtrees across phylogenetic space.

    \b
    For more information, see the homepage at https://github.com/legumeinfo/azulejo

    Written by Joel Berendzen <joelb@ncgr.org>.
    Copyright (C) 2020. National Center for Genome Resources. All rights reserved.
    License: BSD-3-Clause
    """
    options = click_loguru.get_global_options()
    if warnings_as_errors:
        if not options.quiet:
            print(
                "Runtime warnings (e.g., from pandas) will cause exceptions!"
            )
        warnings.filterwarnings("error")
    f"{parallel}"


@cli.command()
@click.option(
    "--force",
    "-f",
    help="Force overwrites of existing binaries.",
    is_flag=True,
    default=False,
)
@click.option(
    "--accept_licenses",
    "-y",
    help="Accept all licenses.",
    is_flag=True,
    default=False,
)
@click.argument("dependencies", nargs=-1)
@click_loguru.init_logger(logfile=False)
def install(dependencies, force, accept_licenses):
    """
    Check for/install binary dependencies.

    \b
    Example:
        azulejo install all

    """
    options = click_loguru.get_global_options()
    installer = DependencyInstaller(
        DEPENDENCY_DICT,
        pkg_name=NAME,
        install_path=INSTALL_PATH,
        force=force,
        accept_licenses=accept_licenses,
        quiet=options.quiet,
    )
    if dependencies == ():
        print(installer.status(exe_paths=True), end="")
    else:
        installer.install_list(dependencies)


# from .analysis import analyze_clusters  # isort:skip
# from .analysis import length_std_dist  # isort:skip
# from .analysis import outlier_length_dist  # isort:skip
# from .analysis import plot_degree_dists  # isort:skip
# from .core import add_singletons  #  isort:skip
# from .core import adjacency_to_graph  #  isort:skip
# from .core import cluster_in_steps  #  isort:skip
# from .core import clusters_to_histograms  #  isort:skip
# from .core import combine_clusters  #  isort:skip
# from .core import compare_clusters  #  isort:skip
# from .core import prepare_protein_files  #  isort:skip
# from .core import homology_cluster  #  isort:skip
from .homology import cluster_build_trees  # isort:skip
from .homology import info_to_fasta  # isort:skip
from .ingest import ingest_sequences  # isort:skip
from .parquet import change_compression  # isort:skip
from .parquet import tsv_to_parquet  # isort:skip
from .proxy import calculate_proxy_genes  # isort:skip
from .synteny import synteny_anchors  # isort:skip
from .synteny import dagchainer_synteny  # isort:skip
from .taxonomy import check_taxonomic_rank  # isort:skip
