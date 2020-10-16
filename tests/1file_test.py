# -*- coding: utf-8 -*-
"""Tests for input table creation."""
# standard library imports
import shutil
from pathlib import Path

# first-party imports
from azulejo.common import dotpath_to_path
from azulejo.ingest import TaxonomicInputTable
from azulejo.ingest import read_from_url

# module imports
from . import W05_INPUTS
from . import W82_INPUTS
from . import print_docstring

# global constants
W05_URL = "https://v1.legumefederation.org/data/index/public/Glycine_soja/W05.gnm1.ann1.T47J/"
SINGLE_INPUT_FILE = "W05.toml"
W82_URL = "https://v1.legumefederation.org/data/index/public/Glycine_max/Wm82.gnm4.ann1.T8TQ/"


@print_docstring()
def test_clean_datadir(request):
    """Clean up datadir."""
    testdir = Path(request.fspath.dirpath())
    datadir = testdir / "data"
    if datadir.exists():
        shutil.rmtree(datadir)  # remove anything left in data directory


@print_docstring()
def test_setup_datadir(request, datadir_mgr, capsys):
    """Copy in and download static data."""
    testdir = Path(request.fspath.dirpath())
    datadir = testdir / "data"
    filesdir = testdir / "testdata"
    shutil.copytree(filesdir, datadir)
    with capsys.disabled():
        datadir_mgr.download(
            download_url=W05_URL,
            files=W05_INPUTS,
            scope="global",
            md5_check=False,
            gunzip=True,
            progressbar=True,
        )
        datadir_mgr.download(
            download_url=W82_URL,
            files=W82_INPUTS,
            scope="global",
            md5_check=False,
            gunzip=True,
            progressbar=True,
        )


@print_docstring()
def test_input_table_parsing(datadir_mgr):
    """Test input table parsing."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=W05_INPUTS + [SINGLE_INPUT_FILE],
        save_outputs=True,
        outscope="global",
        excludepaths=["logs/"],
    ):
        input = TaxonomicInputTable(SINGLE_INPUT_FILE)
        assert input.depth == 3
        assert input.setname == "glycines"
        input_table = input.input_table
        assert len(input_table) == 1
        assert len(input_table.columns) == 7
        rootpath = Path("glycines")
        assert (rootpath / "proteomes.tsv").exists()
        assert (rootpath / "input.toml").exists()
        for subdir in ["glyso", "W05"]:
            rootpath /= subdir
            assert (rootpath / "node_properties.json").exists()
        for i, row in input_table.iterrows():
            out_path = dotpath_to_path(row["path"])
            name = ".".join(row["path"].split(".")[1:])
            for outname, url in [
                (f"{name}.gff", row["gff_url"]),
                (f"{name}.faa", row["fasta_url"]),
            ]:
                with (out_path / outname).open("w") as out_fh:
                    with read_from_url(url) as in_fh:
                        for line in in_fh.read():
                            out_fh.write(line)
