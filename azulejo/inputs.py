# -*- coding: utf-8 -*-
"""Read a TOML taxonomic input table and read files from URLs."""
# standard library imports\
import contextlib
import json
import os
import sys
import shutil
import tempfile
from pathlib import Path

# third-party imports
import attr
import pandas as pd
import smart_open
import toml
from pathvalidate import validate_filename
from pathvalidate import ValidationError
from loguru import logger

# module imports
from .common import dotpath_to_path
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
            input_table_path = root_path / "proteomes.tsv"
            logger.debug(
                f"Input table of {len(self.input_table)} genomes written to {input_table_path}"
            )
            self.input_table.to_csv(input_table_path, sep="\t")
        saved_input_path = root_path / "input.toml"
        if toml_path != saved_input_path:
            shutil.copy2(toml_path, root_path / "input.toml")

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
        # logger.debug(f"Downloading/decompressing {url} to {dirpath}")
        filehandle = smart_open.open(url)
        dldata = filehandle.read()

        def cleanup():
            shutil.rmtree(dirpath)

        with _cd(dirpath, cleanup):

            with open(uncompressed_filename, "w") as f:
                f.write(dldata)
            tmpfile = str(Path(dirpath) / uncompressed_filename)
            yield tmpfile
