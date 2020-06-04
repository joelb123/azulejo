# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# standard library imports
import os
import shutil
from pathlib import Path

# first-party imports
from azulejo.uri import read_input_table

# third-party imports
import pytest
import sh

# global constants
azulejo = sh.Command("azulejo")
DOWNLOAD_URL = "https://v1.legumefederation.org/data/index/public/Glycine_soja/W05.gnm1.ann1.T47J/"
RAW_FASTA_FILE = "glyso.W05.gnm1.ann1.T47J.protein.faa"
RAW_GFF_FILE = "glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3"
INPUT_FILE = "glyso.tsv"


def test_setup(request):
    """Remove datadir, if it exists, and install copies of static data."""
    testdir = Path(request.fspath.dirpath())
    datadir = testdir / "data"
    if datadir.exists():
        shutil.rmtree(datadir)  # remove anything left in data directory
    filesdir = testdir / "testdata"
    shutil.copytree(filesdir, datadir)


def test_input_tsv_parsing(datadir_mgr):
    """Test basic cli function."""
    datadir_mgr.download(
        download_url=DOWNLOAD_URL,
        files=[RAW_FASTA_FILE, RAW_GFF_FILE, INPUT_FILE],
        scope="module",
        md5_check=False,
        gunzip=True,
        progressbar=False,
    )
    with datadir_mgr.in_tmp_dir(
        inpathlist=[RAW_GFF_FILE], save_outputs=True, excludepatterns=["*.log"]
    ):
        outlist = read_input_table(INPUT_FILE)
        print(outlist)
        assert len(outlist) == 2
