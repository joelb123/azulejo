# -*- coding: utf-8 -*-
"""Tests for input table creation."""
# standard library imports
import shutil
from pathlib import Path

# first-party imports
from azulejo.common import dotpath_to_path
from azulejo.ingest import read_from_url
from azulejo.ingest import TaxonomicInputTable


# global constants
DOWNLOAD_URL = "https://v1.legumefederation.org/data/index/public/Glycine_soja/W05.gnm1.ann1.T47J/"
RAW_FASTA_FILE = "glyso.W05.gnm1.ann1.T47J.protein_primaryTranscript.faa"
RAW_GFF_FILE = "glyso.W05.gnm1.ann1.T47J.gene_models_main.gff3"
INPUT_FILE = "glyma+glyso.toml"


def test_setup(request):
    """Remove datadir, if it exists, and install copies of static data."""
    testdir = Path(request.fspath.dirpath())
    datadir = testdir / "data"
    if datadir.exists():
        shutil.rmtree(datadir)  # remove anything left in data directory
    filesdir = testdir / "testdata"
    shutil.copytree(filesdir, datadir)


def test_input_table_downloads(datadir_mgr):
    """Test basic cli function."""
    datadir_mgr.download(
        download_url=DOWNLOAD_URL,
        files=[RAW_FASTA_FILE, RAW_GFF_FILE],
        scope="global",
        md5_check=False,
        gunzip=True,
        progressbar=False,
    )
    with datadir_mgr.in_tmp_dir(
        inpathlist=[RAW_GFF_FILE, RAW_FASTA_FILE, INPUT_FILE],
        save_outputs=True,
        outscope="global",
    ):
        input = TaxonomicInputTable(INPUT_FILE)
        assert input.depth == 3
        assert input.setname == "glycines"
        input_table = input.input_table
        print(f"input_table={input_table}")
        assert len(input_table) == 3
        assert len(input_table.columns) == 7
        rootpath = Path("glycines")
        assert (rootpath / "proteomes.tsv").exists()
        assert (rootpath / "input.toml").exists()
        for subdir in ["glyso", "W05"]:  # sample down depth
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
                        print(f"Copying {url} to {out_path / outname}")
                        for line in in_fh.read():
                            out_fh.write(line)
