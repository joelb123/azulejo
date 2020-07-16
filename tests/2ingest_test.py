# -*- coding: utf-8 -*-
"""Tests for data ingestion."""
# standard library imports
from pathlib import Path

# third-party imports
import sh

# global constants
azulejo = sh.Command("azulejo")
W05_FASTA_FILE = "glyso.W05.gnm1.ann1.T47J.protein.faa"
W05_GFF_FILE = "glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3"
PI_FASTA_FILE = "glycines/glyso/PI483463/glyso.PI483463.faa"
PI_GFF_FILE = "glycines/glyso/PI483463/glyso.PI483463.gff"
LOCAL_INPUT_FILE = "glyso.toml"
NET_INPUT_FILE = "glyma+glyso.toml"


def test_local_data_ingestion(datadir_mgr):
    """Test reading data from local files."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=[
            W05_GFF_FILE,
            W05_FASTA_FILE,
            LOCAL_INPUT_FILE,
            PI_FASTA_FILE,
            PI_GFF_FILE,
        ],
    ):
        output = azulejo(["ingest-sequence-data", LOCAL_INPUT_FILE])
        print(output)
        assert Path("glycines/proteomes.tsv").exists()
        assert Path("glycines/glyso/W05/fragments.tsv").exists()
        assert Path("glycines/glyso/W05/proteins.parq").exists()


def test_net_data_ingestion(datadir_mgr):
    """Test reading compressed data from https."""
    with datadir_mgr.in_tmp_dir(
        inpathlist=[W05_GFF_FILE, W05_FASTA_FILE, NET_INPUT_FILE],
        save_outputs=True,
        outscope="global",
        excludepaths=["logs/"],
    ):
        output = azulejo(["ingest-sequence-data", NET_INPUT_FILE])
        print(output)
        assert Path("glycines/glyma/fragments.tsv").exists()
        assert Path("glycines/glyma/proteins.parq").exists()
