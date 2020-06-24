# -*- coding: utf-8 -*-
"""Tests for basic CLI function."""
# standard library imports
import shutil
from pathlib import Path

# third-party imports
import sh

# global constants
azulejo = sh.Command("azulejo")
W05_FASTA_FILE = "glyso.W05.gnm1.ann1.T47J.protein.faa"
W05_GFF_FILE = "glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3"
PI_FASTA_FILE = "glycines/glyso/PI483463/glyso.PI483463.faa"
PI_GFF_FILE = "glycines/glyso/PI483463/glyso.PI483463.gff"
INPUT_FILE = "glyso.toml"


def test_local_data_ingestion(datadir_mgr):
    """Test basic cli function."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=[
            W05_GFF_FILE,
            W05_FASTA_FILE,
            INPUT_FILE,
            PI_FASTA_FILE,
            PI_GFF_FILE,
        ],
        save_outputs=True,
        outscope="global",
    ):
        output = azulejo(["ingest-data", "glyso.toml"])
        print(output)
