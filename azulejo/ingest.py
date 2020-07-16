# -*- coding: utf-8 -*-
"""Sequence (FASTA) and genome (GFF) ingestion operations."""
# standard library imports
import contextlib
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

# third-party imports
import attr
import click
import dask.bag as db
import gffpandas.gffpandas as gffpd
import pandas as pd
import toml
from dask.diagnostics import ProgressBar

# first-party imports
import smart_open
from loguru import logger
from pathvalidate import ValidationError
from pathvalidate import validate_filename

# module imports
from . import cli
from . import click_loguru
from .common import FRAGMENTS_FILE
from .common import PROTEINS_FILE
from .common import PROTEOMES_FILE
from .common import SAVED_INPUT_FILE
from .common import dotpath_to_path
from .common import sort_proteome_frame
from .common import write_tsv_or_parquet
from .core import cleanup_fasta
from .taxonomy import rankname_to_number

# global constants
__all__ = ["TaxonomicInputTable", "read_from_url"]
FILE_TRANSPORT = "file://"
REQUIRED_LEAF_NAMES = (
    "fasta",
    "gff",
)
COMPRESSION_EXTENSIONS = (
    "gz",
    "bz2",
)


@cli.command()
@click_loguru.init_logger()
@click.option(
    "-s/-l",
    "--shorten_source/--no-shorten_source",
    default=True,
    is_flag=True,
    show_default=True,
    help="Remove invariant dotpaths in source IDs.",
)
@click.option(
    "--parallel/--no-parallel",
    is_flag=True,
    default=True,
    show_default=True,
    help="Process in parallel.",
)
@click.argument("input_toml")
def ingest_sequence_data(shorten_source, input_toml, parallel):
    """Marshal protein and genome sequence information.

    IDs must correspond between GFF and FASTA files and must be unique across
    the entire set.
    """
    options = click_loguru.get_global_options()
    input_obj = TaxonomicInputTable(Path(input_toml), write_table=False)
    input_table = input_obj.input_table
    set_path = Path(input_obj.setname)
    arg_list = []
    for i, row in input_table.iterrows():
        arg_list.append((row["path"], row["fasta_url"], row["gff_url"],))
    if parallel:
        bag = db.from_sequence(arg_list)
    else:
        file_stats = []
    if not options.quiet:
        logger.info(f"Extracting FASTA/GFF info for {len(arg_list)} genomes:")
        ProgressBar().register()
    if parallel:
        file_stats = bag.map(
            read_fasta_and_gff, shorten_source=shorten_source
        ).compute()
    else:
        for args in arg_list:
            file_stats.append(
                read_fasta_and_gff(args, shorten_source=shorten_source)
            )
    del arg_list
    seq_frame = pd.DataFrame.from_dict([s[0] for s in file_stats]).set_index(
        "path"
    )
    frag_frame = pd.DataFrame.from_dict([s[1] for s in file_stats]).set_index(
        "path"
    )
    proteomes = pd.concat(
        [input_table.set_index("path"), frag_frame, seq_frame], axis=1
    )
    proteomes.drop(["fasta_url", "gff_url"], axis=1, inplace=True)
    proteomes = sort_proteome_frame(proteomes)
    proteome_table_path = set_path / PROTEOMES_FILE
    logger.info(
        f'Writing table of proteomes to "{proteome_table_path}", edit it to'
        " change preferences"
    )
    logger.info(
        "This is also the time to put DNA fragments on a common naming basis."
    )
    write_tsv_or_parquet(proteomes, proteome_table_path)


def read_fasta_and_gff(args, shorten_source=True):
    """Read corresponding sequence and position files and construct consolidated tables."""
    dotpath, fasta_url, gff_url = args
    out_path = dotpath_to_path(dotpath)

    with read_from_url(fasta_url) as fasta_fh:
        unused_stem, unused_path, prop_frame, file_stats = cleanup_fasta(
            out_path, fasta_fh, dotpath, write_fasta=False, write_stats=False
        )
    # logger.debug(f"Reading GFF file {gff_url}.")
    with filepath_from_url(gff_url) as local_gff_file:
        annotation = gffpd.read_gff3(local_gff_file)
    features = annotation.filter_feature_of_type(
        ["mRNA"]
    ).attributes_to_columns()
    del annotation
    features.drop(
        features.columns.drop(
            ["seq_id", "start", "strand", "ID"]
        ),  # drop EXCEPT these
        axis=1,
        inplace=True,
    )  # drop non-essential columns
    if shorten_source:
        # drop identical sub-fields in seq_id to make them easier to map
        split_sources = features["seq_id"].str.split(".", expand=True)
        split_sources = split_sources.drop(
            [
                i
                for i in split_sources.columns
                if len(set(split_sources[i])) == 1
            ],
            axis=1,
        )
        sources = split_sources.agg(".".join, axis=1)
        features["seq_id"] = sources
    features = features.set_index("ID")
    # Make a categorical column, frag_id, based on seq_id
    features["frag_id"] = pd.Categorical(features["seq_id"])
    features.drop(["seq_id"], axis=1, inplace=True)
    # Drop any features not found in sequence file, e.g., zero-length
    features = features[features.index.isin(prop_frame.index)]
    # sort fragments by largest value
    frag_counts = features["frag_id"].value_counts()
    frag_frame = pd.DataFrame()
    frag_frame["counts"] = frag_counts
    frag_frame["idx"] = range(len(frag_frame))
    frag_frame["frag_id"] = frag_frame.index
    frag_frame["new_name"] = ""
    frag_frame.set_index(["idx"], inplace=True)
    frag_stats = {
        "path": dotpath,
        "frag.n": len(frag_frame),
        "frag.max": frag_counts[0],
    }
    frag_count_path = out_path / FRAGMENTS_FILE
    if not frag_count_path.exists():
        # logger.debug(f"Writing fragment stats to file {frag_count_path}")
        write_tsv_or_parquet(frag_frame, frag_count_path)
    del frag_frame
    features["frag_count"] = features["frag_id"].map(frag_counts)
    features.sort_values(
        by=["frag_count", "start"], ascending=[False, True], inplace=True
    )
    frag_id_range = []
    for frag_id in frag_counts.index:
        frag_id_range += list(range(frag_counts[frag_id]))
    features["frag_pos"] = frag_id_range
    del frag_id_range
    features.drop(["frag_count"], axis=1, inplace=True)
    # join GFF info to FASTA info
    joined_path = out_path / PROTEINS_FILE
    features = features.join(prop_frame)
    write_tsv_or_parquet(features, joined_path)
    return file_stats, frag_stats


class TaxonomicInputTable:

    """Parse an azulejo input dictionary."""

    def __init__(self, toml_path, write_table=True):
        """Create structures."""
        self.depth = 0
        try:
            tree = toml.load(toml_path)
        except TypeError:
            logger.error(f"Error in filename {toml_path}")
            sys.exit(1)
        except toml.TomlDecodeError as e:
            logger.error(f"File {toml_path} is not valid TOML")
            logger.error(e)
            sys.exit(1)
        if len(tree) > 1:
            logger.error(
                f"Input file {toml_path} should define a single "
                + f"object, but defines {len(tree)} instead"
            )
            sys.exit(1)
        self.setname = self._validate_name(list(tree.keys())[0])
        root_path = Path(self.setname)
        if not root_path.exists():
            logger.info(f"Creating directory for set {self.setname}")
            root_path.mkdir(parents=True)
        self._Node = attr.make_class(
            "Node", ["path", "name", "rank", "rank_val"]
        )
        self._nodes = []
        self._genome_dir_dict = {}
        self._n_genomes = 0
        self._walk(self.setname, tree[self.setname])
        self.input_table = pd.DataFrame.from_dict(
            self._genome_dir_dict
        ).transpose()
        del self._genome_dir_dict
        del self._nodes
        self.input_table.index.name = "order"
        if write_table:
            input_table_path = root_path / PROTEOMES_FILE
            logger.debug(
                f"Input table of {len(self.input_table)} genomes written to"
                f" {input_table_path}"
            )
            write_tsv_or_parquet(self.input_table, input_table_path)
        saved_input_path = root_path / SAVED_INPUT_FILE
        if toml_path != saved_input_path:
            shutil.copy2(toml_path, root_path / SAVED_INPUT_FILE)

    def _validate_name(self, name):
        """Check if a potential filename is valid or not."""
        try:
            validate_filename(name)
        except ValidationError as e:
            logger.error(f"Invalid component name {name} in input file")
            sys.exit(1)
        return name

    def _validate_uri(self, uri):
        """Check if the transport at the start of a URI is valid or not."""
        try:
            smart_open.parse_uri(uri)
        except NotImplementedError:
            logger.error(f'Unimplemented transport in uri "{uri}"')
            sys.exit(1)
        return uri

    def _strip_file_uri(self, url):
        """Removes the file:// uri from a URL string."""
        if url.startswith(FILE_TRANSPORT):
            return url[len(FILE_TRANSPORT) :]
        return url

    def _walk(self, node_name, tree):
        """Recursively walk tree structure."""
        # Check for required field properties.
        if len(self._nodes) > 0:
            dot_path = f"{self._nodes[-1].path}.{node_name}"
        else:
            dot_path = node_name
        if "name" not in tree:
            tree["name"] = f"'{node_name}'"
        if "rank" not in tree:
            logger.error(f'Required entry "rank" not in entry {dot_path}')
            sys.exit(1)
        try:
            rank_val = rankname_to_number(tree["rank"])
        except ValueError as e:
            logger.error(f"Unrecognized taxonomic rank {tree['rank']}")
            logger.error(e)
            sys.exit(1)
        if (len(self._nodes) > 0) and rank_val <= self._nodes[-1].rank_val:
            logger.error(
                f"rank {tree['rank']} value {rank_val} is not less than"
                + f" previous rank value of {self._nodes[-1].rank_val}"
            )
            sys.exit(1)
        # Push node onto stack
        self._nodes.append(
            self._Node(
                self._validate_name(dot_path),
                tree["name"],
                tree["rank"],
                rank_val,
            )
        )
        self.depth = max(self.depth, len(self._nodes))
        # Initialize node properties dictionary
        properties = {"path": dot_path, "children": []}
        for k, v in tree.items():
            if isinstance(v, dict):
                properties["children"] += [k]
                self._walk(k, v)
            else:
                properties[k] = v
        if len(properties["children"]) == 0:
            del properties["children"]
        # Check if this is a genome directory node
        genome_dir_properties = [
            (p in properties) for p in REQUIRED_LEAF_NAMES
        ]
        if any(genome_dir_properties):
            if not all(genome_dir_properties):
                missing_properties = [
                    p
                    for i, p in enumerate(REQUIRED_LEAF_NAMES)
                    if not genome_dir_properties[i]
                ]
                logger.error(
                    f"Missing properties {missing_properties} "
                    + f"for node {dot_path}"
                )
                sys.exit(1)
            if "uri" not in tree:
                uri = FILE_TRANSPORT
            else:
                uri = self._validate_uri(tree["uri"])
                if not uri.endswith("/"):
                    uri += "/"
            self._genome_dir_dict[self._n_genomes] = {"path": dot_path}
            if "preference" not in tree:
                self._genome_dir_dict[self._n_genomes]["preference"] = ""
            else:
                self._genome_dir_dict[self._n_genomes]["preference"] = tree[
                    "preference"
                ]
            for n in self._nodes:
                self._genome_dir_dict[self._n_genomes][
                    f"phy.{n.rank}"
                ] = n.name
            self._genome_dir_dict[self._n_genomes][
                "fasta_url"
            ] = self._strip_file_uri(uri + tree["fasta"])
            self._genome_dir_dict[self._n_genomes][
                "gff_url"
            ] = self._strip_file_uri(uri + tree["gff"])
            self._n_genomes += 1
        for n in self._nodes:
            properties[n.rank] = n.name
        node_path = dotpath_to_path(dot_path)
        node_path.mkdir(parents=True, exist_ok=True)
        properties_file = node_path / "node_properties.json"
        logger.debug(f"Writing properties file to {properties_file}")
        with properties_file.open("w") as filepointer:
            json.dump(properties, filepointer)
        # Pop node from stack
        self._nodes.pop()


@contextlib.contextmanager
def _cd(newdir, cleanup=lambda: True):
    "Change directory with cleanup."
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)
        cleanup()


@contextlib.contextmanager
def read_from_url(url):
    """Read from a URL transparently decompressing if compressed."""
    yield smart_open.open(url)


@contextlib.contextmanager
def filepath_from_url(url):
    """Get a local file from a URL, decompressing if needed."""
    filename = url.split("/")[-1]
    compressed = False
    uncompressed_filename = filename
    for ext in COMPRESSION_EXTENSIONS:
        if filename.endswith(ext):
            compressed = True
            uncompressed_filename = filename[: -(len(ext) + 1)]
            break
    if (
        url.find("://") == -1 and not compressed
    ):  # no transport, must be a file
        yield url
    else:
        dirpath = tempfile.mkdtemp()
        filehandle = smart_open.open(url)
        dldata = filehandle.read()

        def cleanup():
            shutil.rmtree(dirpath)

        with _cd(dirpath, cleanup):

            with open(uncompressed_filename, "w") as f:
                f.write(dldata)
            tmpfile = str(Path(dirpath) / uncompressed_filename)
            yield tmpfile
