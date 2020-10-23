# -*- coding: utf-8 -*-
"""azulejo -- tile phylogenetic space with subtrees."""
# standard library imports
import locale
import os
import sys
import warnings
from pathlib import Path
from pkg_resources import iter_entry_points

# third-party imports
import click
from click_plugins import with_plugins
from click_loguru import ClickLoguru

# module imports
from .common import NAME
from .core import homology_cluster as undeco_homology_cluster
from .core import cluster_in_steps as undeco_cluster_in_steps
from .homology import cluster_build_trees as undeco_cluster_build_trees
from .homology import info_to_fasta as undeco_into_to_fasta
from .ingest import ingest_sequences as undeco_ingest_sequences
from .installer import DependencyInstaller
from .parquet import parquet_to_tsv as undeco_parquet_to_tsv
from .proxy import calculate_proxy_genes as undeco_calculate_proxy_genes
from .synteny import intersect_anchors as undeco_intersect_anchors
from .synteny import synteny_anchors as undeco_synteny_anchors
from .synteny import unique_anchors as undeco_unique_anchors
from .taxonomy import print_taxonomic_ranks
from .taxonomy import rankname_to_number

# global constants
LOG_FILE_RETENTION = 3
__version__ = "0.9.14"
INSTALL_ENVIRON_VAR = (  # installs go into "/bin" and other subdirs of this directory
    NAME.upper() + "_INSTALL_DIR"
)
if INSTALL_ENVIRON_VAR in os.environ:
    INSTALL_PATH = Path(os.environ[INSTALL_ENVIRON_VAR])
else:
    INSTALL_PATH = None
MUSCLE_VER = "3.8.1551"
USEARCH_VER = "11.0.667"
DAGCHAINER_TOOL_VER = "1.0"
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
        "version_parser": lambda v: v.split()[1],
        "make": ["muscle-pgo"],
        "copy_binaries": ["muscle"],
        "license": "public domain",
        "license_restrictive": False,
    },
    "usearch": {
        "binaries": ["usearch"],
        "dir": ".",
        "version": USEARCH_VER,
        "version_command": ["-version"],
        "version_parser": lambda v: v.split()[1].split("_")[0],
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
    "dagchainer_tool": {
        "binaries": ["dagchainer_tool.sh"],
        "tarball": (
            "https://sourceforge.net/projects/dagchainer/files/dagchainer/"
            + "DAGchainer-r02062008/DAGchainer_r02-06-2008.tar.gz/download"
        ),
        "dir": "DAGCHAINER",
        "required": False,
        "version": DAGCHAINER_TOOL_VER,
        "version_command": ["version"],
        "version": "1.0",
        "patch": ["-p1", "-i", "dagchainer-0206008.patch"],
        "make": [],
        "copy_binaries": [
            "dagchainer",
            "dagchainer_tool.sh",
            "hash_into_fasta_id.pl",
            "pairs_to_adjacency.py",
            "top_blast_hit.awk",
            "run_DAG_chainer.pl",
        ],
        "license": "GPL-2",
        "license_restrictive": False,
    },
    "blast": {
        "binaries": ["blastx"],
        "tarball": (
            "https://sourceforge.net/projects/dagchainer/files/dagchainer/"
            + "DAGchainer-r02062008/DAGchainer_r02-06-2008.tar.gz/download"
        ),
        "dir": "DAGCHAINER",
        "required": False,
        "version": DAGCHAINER_TOOL_VER,
        "version_command": ["-version"],
        "version": "2.10.1",
        "version_parser": lambda v: v.split()[1].split("_")[0].rstrip("+"),
        "patch": ["-p1", "-i", "dagchainer-0206008.patch"],
        "make": [],
        "copy_binaries": [
            "dagchainer",
            "dagchainer_tool.sh",
            "blinkPerl_v1.1.pl",
            "hash_into_fasta_id.pl",
            "pairs_to_adjacency.py",
            "rename_reorder_adjacency.py",
            "top_blast_hit.awk",
        ],
        "license": "GPL-2",
        "license_restrictive": False,
    },
}
DEFAULT_K = 2
DEFAULT_STEPS = 16

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
def cli(warnings_as_errors, parallel, **kwargs):
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
    f"{parallel}{kwargs}"


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
    dep_installer = DependencyInstaller(
        DEPENDENCY_DICT,
        pkg_name=NAME,
        install_path=INSTALL_PATH,
        force=force,
        accept_licenses=accept_licenses,
        quiet=options.quiet,
    )
    if dependencies == ():
        print(dep_installer.status(exe_paths=True), end="")
    else:
        dep_installer.install_list(dependencies)


@cli.command()
@click_loguru.init_logger()
@click_loguru.log_elapsed_time(level="info")
@click_loguru.log_peak_memory_use(level="info")
@click.option(
    "--identity",
    "-i",
    default=0.0,
    help="Minimum sequence ID (0-1). [default: lowest]",
)
@click.option(
    "--cluster_file",
    "-c",
    type=click.Path(exists=True),
    default=None,
    show_default=True,
    help="Use pre-existing homology clusters.",
)
@click.argument("setname")
def homology(identity, setname, cluster_file):
    """
    Calculate homology clusters, MSAs, trees.

    \b
    Example:
        azulejo homology glycines

    """
    undeco_cluster_build_trees(
        identity, setname, cluster_file=cluster_file, click_loguru=click_loguru
    )


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--append/--no-append",
    "-a/-x",
    is_flag=True,
    default=True,
    help="Append to FASTA file.",
    show_default=True,
)
@click.argument("parquetfile")
@click.argument("fastafile")
def parquet_to_fasta(parquetfile, fastafile, append):
    """Convert Parquet sequence info to FASTA file."""
    undeco_into_to_fasta(parquetfile, fastafile, append)


@cli.command()
@click_loguru.init_logger()
@click_loguru.log_elapsed_time(level="info")
@click_loguru.log_peak_memory_use(level="info")
@click.option(
    "-k", default=DEFAULT_K, help="Synteny anchor length.", show_default=True
)
@click.option(
    "--peatmer/--kmer",
    default=True,
    is_flag=True,
    show_default=True,
    help="Allow repeats in anchor.",
)
@click.option(
    "--write_ambiguous/--no-write-ambiguous",
    default=True,
    is_flag=True,
    show_default=True,
    help="Include ambiguous anchors.",
)
@click.argument("setname")
def synteny(k, peatmer, setname, write_ambiguous):
    """Calculate synteny anchors.

    \b
    Example:
        azulejo synteny glycines

    """
    undeco_synteny_anchors(
        k,
        peatmer,
        setname,
        click_loguru=click_loguru,
        write_ambiguous=write_ambiguous,
    )


@cli.command()
@click_loguru.init_logger()
@click_loguru.log_elapsed_time(level="info")
@click.option(
    "-k", default=DEFAULT_K, help="Synteny anchor length.", show_default=True
)
@click.argument("setname")
def unique_anchors(setname, k):
    """Uniqueify synteny anchors.

    \b
    Example:
        azulejo unique_anchors glycines

    """
    undeco_unique_anchors(
        setname, k,
    )


@cli.command()
@click_loguru.init_logger()
@click_loguru.log_elapsed_time(level="info")
@click.option(
    "-k", default=DEFAULT_K, help="Synteny anchor length.", show_default=True
)
@click.argument("setname")
@click.argument("compfile")
def intersect_anchors(setname, k, compfile):
    """Intersect two sets of synteny anchors.

    \b
    Example:
        azulejo intersect_anchors glycines dagchainer_tool_out/synteny_anchor_summary.tsv

    """
    undeco_intersect_anchors(setname, k, compfile)


@cli.command()
@click.argument("rankname", nargs=-1)
def taxonomy(rankname):
    """Check/show taxonomic ranks."""
    if rankname == ():
        print_taxonomic_ranks()
    else:
        rankname = rankname[0]
        try:
            rankval = rankname_to_number(rankname)
        except ValueError as error_msg:
            print(f"ERROR: {error_msg}")
            sys.exit(1)
        print(rankval)


@cli.command()
@click_loguru.init_logger(logfile=False)
@click.option(
    "--columns",
    default=False,
    is_flag=True,
    show_default=False,
    help="Print names/dtypes of columns.",
)
@click.option(
    "--index_name",
    default=False,
    is_flag=True,
    show_default=False,
    help="Print the name of the index.",
)
@click.option(
    "--no_index",
    default=False,
    is_flag=True,
    show_default=False,
    help="Do not write index column.",
)
@click.option(
    "--pretty",
    default=False,
    is_flag=True,
    show_default=False,
    help="Pretty-print output.",
)
@click.option(
    "--no_header",
    "-h",
    default=False,
    is_flag=True,
    show_default=False,
    help="Do not write the header.",
)
@click.option(
    "--col",
    "-c",
    default=None,
    multiple=True,
    show_default=True,
    help="Write only the named column.",
)
@click.option(
    "--writefile",
    "-w",
    default=False,
    is_flag=True,
    show_default=True,
    help="Write to a TSV file.",
)
@click.option(
    "--index_val",
    "-i",
    default=None,
    multiple=True,
    show_default=True,
    help="Write only the row with this index value.",
)
@click.option(
    "--head",
    default=0,
    show_default=False,
    help="Write only the first N rows.",
)
@click.option(
    "--max_rows",
    default=None,
    show_default=False,
    help="Pretty-print N rows.",
)
@click.option(
    "--max_cols",
    default=None,
    show_default=False,
    help="Pretty-print N cols.",
)
@click.option(
    "--tail",
    default=0,
    show_default=False,
    help="Write only the last N rows.",
)
@click.argument("parquetfile", type=click.Path(exists=True))
@click.argument("tsvfile", nargs=-1)
def parquet_to_tsv(
    parquetfile,
    tsvfile,
    columns,
    index_name,
    col,
    no_index,
    head,
    tail,
    index_val,
    no_header,
    pretty,
    max_rows,
    max_cols,
    writefile,
):
    """
    Reads parquet file, writes tsv.

    \b
    Example:
        azulejo parquet-to-tsv glycines/proteomes.syn.parq
    """
    undeco_parquet_to_tsv(
        parquetfile,
        tsvfile,
        columns,
        index_name,
        col,
        no_index,
        head,
        tail,
        index_val,
        no_header,
        pretty,
        max_rows,
        max_cols,
        writefile,
    )


@cli.command()
@click_loguru.init_logger()
@click.argument("seqfile")
@click.option(
    "--identity",
    "-i",
    default=0.0,
    help="Minimum sequence identity (float, 0-1). [default: lowest]",
)
@click.option(
    "--min_id_freq",
    "-m",
    default=0,
    show_default=True,
    help="Minimum frequency of ID components.",
)
@click.option(
    "--delete/--no-delete",
    "-d/-n",
    is_flag=True,
    default=True,
    help="Delete primary output of clustering. [default: delete]",
)
@click.option(
    "--write_ids/--no-write_ids",
    "-w",
    is_flag=True,
    default=False,
    help="Write file of ID-to-clusters. [default: delete]",
)
@click.option(
    "--do_calc/--no-do_calc",
    "-c/-x",
    is_flag=True,
    default=True,
    help="Write file of ID-to-clusters. [default: do_calc]",
)
@click.option(
    "--substrs", help="subpath to file of substrings. [default: none]"
)
@click.option("--dups", help="subpath to file of duplicates. [default: none]")
def cluster(
    seqfile,
    identity,
    delete=True,
    write_ids=False,
    do_calc=True,
    min_id_freq=0,
    substrs=None,
    dups=None,
    cluster_stats=True,
    outname=None,
):
    """Cluster at a global sequence identity threshold."""
    undeco_homology_cluster(
        seqfile,
        identity,
        delete=delete,
        write_ids=write_ids,
        do_calc=do_calc,
        min_id_freq=min_id_freq,
        substrs=substrs,
        dups=dups,
        cluster_stats=cluster_stats,
        outname=outname,
        click_loguru=click_loguru,
    )


@cli.command()
@click_loguru.init_logger()
@click.argument("seqfile")
@click.option(
    "--steps",
    "-s",
    default=DEFAULT_STEPS,
    show_default=True,
    help="# of steps from lowest to highest",
)
@click.option(
    "--min_id_freq",
    "-m",
    default=0,
    show_default=True,
    help="Minimum frequency of ID components.",
)
@click.option(
    "--substrs", help="subpath to file of substrings. [default: none]"
)
@click.option("--dups", help="subpath to file of duplicates. [default: none]")
def cluster_in_steps(seqfile, steps, min_id_freq=0, substrs=None, dups=None):
    """Cluster in steps from low to 100% identity."""
    undeco_cluster_in_steps(
        seqfile, steps, min_id_freq=min_id_freq, substrs=substrs, dups=dups
    )


@cli.command()
@click_loguru.init_logger()
@click_loguru.log_elapsed_time()
@click.argument("input_toml")
def ingest(input_toml):
    """
    Marshal protein and genome sequence information.

    IDs must correspond between GFF and FASTA files and must be unique across
    the entire set.

    \b
    Example:
        azulejo ingest glyma+glyso.toml

    """
    undeco_ingest_sequences(input_toml, click_loguru=click_loguru)


@cli.command()
@click_loguru.init_logger()
@click.argument("setname")
@click.argument("synteny_type")
@click.argument("prefs", nargs=-1)
def proxy_genes(setname, synteny_type, prefs):
    """Calculate a set of proxy genes."""
    undeco_calculate_proxy_genes(setname, synteny_type, prefs)
