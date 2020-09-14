.. epigraph:: azulejo
              noun INFORMAL
              a glazed tile, usually blue, found on the inside of churches and palaces in Spain and Portugal.

azulejo
=======
``azulejo`` azulejo tiles phylogenetic space with subtrees
normalizes and validates genomic data files prior to further processing
or inclusion in a data store such as that of the
`Legume Federation <https://www.legumefederation.org/en/data-store/>`_.

Prerequisites
-------------
Python 3.6 or greater is required. This package is tested under Linux using Python 3.8.

``azulejo`` relies on `usearch <https://www.drive5.com/usearch/download.html>`_ for 
clustering and `MUSCLE <https://www.drive5.com/muscle/downloads.htm>`_ for sequence
alignment and initial tree-building.  ``usearch`` is free for individual, non-commercial
use, while ``MUSCLE`` is in the public domain.  Both can be downloaded and installed
by the ``azulejo install`` subcommand.  You will be asked to accept the license terms
for ``usearch`` by the install procedure.

Installation for Users
----------------------
Install via pip or (better yet) `pipx <https://pipxproject.github.io/pipx/>`_: ::

     pipx install azulejo

``azulejo`` contains some long commands and many options.  To enable command-line
completion for ``azulejo`` commands, execute the following command if you are using
``bash`` as your shell: ::

    eval "$(_AZULEJO_COMPLETE=source_bash azulejo)"

Then you should run ``azulejo install`` to see if you have the binary dependencies
installed.  If not, ``azulejo install all`` will install them for you.

For Developers
--------------
If you plan to develop ``azulejo``, you'll need to install
the `poetry <https://python-poetry.org>`_ dependency manager.
If you haven't previously installed ``poetry``, execute the command: ::

    curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python

Next, get the master branch from GitHub ::

	git clone https://github.com/legumeinfo/azulejo.git

Change to the ``azulejo/`` directory and install with poetry: ::

	poetry install -v

Run ``azulejo`` with ``poetry``: ::

    poetry run azulejo

Usage
-----
Installation puts a single script called ``azulejo`` in your path.  The usage format is::

    azulejo [GLOBALOPTIONS] COMMAND [COMMANDOPTIONS][ARGS]


Master Input File
-----------------
``azulejo`` uses a configuration file in `TOML  <https://github.com/toml-lang/toml>`_
format as the master input that associates files with phylogeny.  The format of this file
is the familiar headings in square brackets followed by configuration values::

    [glycines]
    rank = "genus"
    name = "Glycine"

    [glycines.glyso]
    rank = "species"
    name = "Glycine soja"

    [glycines.glyso.PI483463]
    rank = "strain"
    gff = "glyso.PI483463.gnm1.ann1.3Q3Q.gene_models_main.gff3.gz"
    fasta = "glyso.PI483463.gnm1.ann1.3Q3Q.protein_primaryTranscript.faa.gz"
    uri = "https://v1.legumefederation.org/data/index/public/Glycine_soja/PI483463.gnm1.ann1.3Q3Q/"
    comments = """
    Glycine soja accession PI 483463 has been identified as being unusually
    salt-tolerant (Lee et al., 2009)."""


* [headings]
    There can be only one top-level heading, and that will be the name of the
    resulting output set.  This name will be the name of an output directory that will be
    created in the current working directory, so this heading (and all subheadings) must
    obey UNIX filesystem naming rules or an error will result.  Each heading level
    (indicated by a ".") will result in another taxonomic level and another directory level
    in the output directory.  Depths do not need to be consistent.

* rank
    Each level defined must have a ``rank`` defined, and that rank must match one of the
    taxonomic ranks defined by ``azulejo``, which you can view and test using the
    ``check-taxonomic-rank`` command.   There are 24 major taxonomic ranks, each of which
    may be modified by 16 different prefixes for a total of 174 taxonomic levels (some of
    which are synonoymous).

* name
    Each level may (and usually should) have a ``name`` defined.  This name is intended
    to be human-readable with no restrictions on the characters used, but it goes into
    plot legends in places, so it's best to not make it too long. If the name is not specified,
    it will be taken from the level name enclosed in single quotes (e.g., 'PI483463' for the
    example above).

* fasta
    If the level specifies a genome, it must have a ``fasta`` entry corresponding
    to the name of the *protein* FASTA file.  In eukaryotes, the FASTA file should be a
    file of primary (generally longest) protein transcripts, if available, rather than all protein
    transcripts (i.e., not including splice variants). Sequences will be cleaned of dashes, stops,
    and other out-of-alphabet characters.  Ambiguous residues at the beginnings and ends of
    sequences will be trimmed. Zero-length sequences will be discarded, which can result in a
    smaller number of sequences out.  These files may be compressed, with extensions ``.gz`` or
    ``.bz2``.

* gff
    If the level specifies a genome, it must have a ``gff`` entry corresponding
    to a version 3 Genome Feature File (GFF3) containing ``CDS`` entries with ID values
    matching those IDs in the FASTA file.  The same compression extensions as for
    ``fasta`` entries apply.  If the ``SOURCE`` fields in those CDS entries
    (which contain the names of the DNA fragments such as scaffolds that the CDS came from)
    contain dot-separated components, those components that are identical across the entire
    file will be discarded by default.  There is an opportunity later in the process to
    remap DNA source names to a common dictionary for comparison among chromosomes and
    plastids.

* uri
    This optional field may contain a a uniform resource identifier such as
    ``https://sitename/dir/``.  ``azulejo`` uses `smart-open <https://www.pypi.org/project/smart-open/>`_
    for doing transparent on-the-fly decompression from a variety of file systems
    including HTTPS, HDFS, SSH, and SFTP (but not FTP).
    If this field is not supplied, local file access is assumed with paths relative to
    the current working directory. The URI will be prepended to ``fasta``
    and ``gff`` paths, allowing for convenient downloading on-the-fly from sites such as
    LegumeInfo or GenBank.   Downloads are not cached, so if you intend to run ``azulejo``
    multiple times on the same input data, you will save time by downloading and uncompressing
    files to local storage.

* preference
    This optional field may be used to override the genome preference heuristic
    that is the fall-thru preference after proxy-gene heuristics have been applied.  This is an integer
    value, with lower integers getting the highest priority.  Set this value to zero if you
    know in advance that one of the input genomes is considered the reference genome and,
    all things being equal, you would prefer to select proxy genes from this genome.  You
    may also set these preference values later, after the default genome preference (genomes
    will be preferred in order of the most genes in a single DNA fragment) has already been
    applied, but before proxy gene selection.

* other info
    A design goal for ``azulejo`` was to not lose metadata, even if it
    was not used by ``azulejo`` itself, while keeping metadata out of file names.
    As an aid in that goal, for each (sub)heading level/output directory, ``azulejo``
    creates a JSON file named ``node_properties.json`` at each node in the output
    hierarchy that containing all information from this file as well as other information
    calculated at ingestion time by ``azulejo``.  You may specify any additional data you would
    like to pass along (e.g., for later use in a web page) and it will be translated from TOML
    to JSON and passed along, such as the multi-line ``comments`` field in the example.
    Examples of useful metadata that may be easier to enter at ingestion time than to
    garner later include taxon IDs of the level and its parent, common names, URLs of
    papers describing the genome, and geographic origin of the sample.

A copy of the input file will be saved in the output directory under the name ``input.toml``.
See the examples in the ``tests/testdata`` repository directory for examples of input data.

Global Options
--------------
The following options are global in scope and, if used must be placed before
``COMMAND``:

============================= ===========================================
   -v, --verbose              Log debugging info to stderr.
   -q, --quiet                Suppress logging to stderr.
   --no-logfile               Suppress logging to file.
   -e, --warnings_as_errors   Treat warnings as fatal (for testing).
============================= ===========================================

Commands
--------
A listing of commands is available via ``azulejo --help``.
The currently implemented commands are, in the order they will normally be run:

========================= ==================================================
  install                 Check for/install binary dependencies.
  ingest                  Marshal protein and genome sequence information.
  homology                Calculate homology clusters, MSAs, trees.
  synteny                 Calculate synteny anchors.
  proxy-genes             Calculate a set of proxy genes from synteny files.
  parquet-to-tsv          Reads parquet file, writes tsv.
========================= ==================================================

``azulejo`` stores most intermediate results in the Parquet format with
extension ``.parq``.  These binary files are compressed and typically can
be read more than 30X faster than the tab-separated-value (TSV) files they
can be interconverted with.  In addition, Parquet files do not lose metadata
such as binary representation sizes.

Each command has its ``COMMANDOPTIONS``, which may be listed with: ::

    azulejo COMMAND --help

Project Status
--------------
+-------------------+-------------+------------+
| Latest Release    | |pypi|      | |azulejo|  |
+-------------------+-------------+            +
| Activity          | |repo|      |            |
+-------------------+-------------+            +
| Downloads         | |downloads| |            |
+-------------------+-------------+            +
| Download Rate     | |dlrate|    |            |
+-------------------+-------------+            +
| License           | |license|   |            |
+-------------------+-------------+            +
| Code Grade        | |codacy|    |            |
+-------------------+-------------+            +
| Coverage          | |coverage|  |            |
+-------------------+-------------+            +
| Travis Build      | |travis|    |            |
+-------------------+-------------+            +
| Issues            | |issues|    |            |
+-------------------+-------------+            +
| Code Style        | |black|     |            |
+-------------------+-------------+------------+


.. |azulejo| image:: docs/azulejo.jpg
     :target: https://en.wikipedia.org/wiki/Azulejo
     :alt: azulejo Definition

.. |black| image:: https://img.shields.io/badge/code%20style-black-000000.svg?style=flat-square
    :target: https://github.com/psf/black
    :alt: Black is the uncompromising Python code formatter

.. |pypi| image:: https://img.shields.io/pypi/v/azulejo.svg
    :target: https://pypi.python.org/pypi/azulejo
    :alt: Python package

.. |repo| image:: https://img.shields.io/github/last-commit/legumeinfo/azulejo
    :target: https://github.com/legumeinfo/azulejo
    :alt: GitHub repository

.. |license| image:: https://img.shields.io/badge/License-BSD%203--Clause-blue.svg
    :target: https://github.com/legumeinfo/azulejo/blob/master/LICENSE
    :alt: License terms

.. |rtd| image:: https://readthedocs.org/projects/azulejo/badge/?version=latest
    :target: http://azulejo.readthedocs.io/en/latest/?badge=latest
    :alt: Documentation Server

.. |travis| image:: https://img.shields.io/travis/legumeinfo/azulejo.svg
    :target:  https://travis-ci.org/legumeinfo/azulejo
    :alt: Travis CI

.. |codacy| image:: https://api.codacy.com/project/badge/Grade/99549f0ed4e6409e9f5e80a2c4bd806b
    :target: https://www.codacy.com/app/joelb123/azulejo?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=legumeinfo/azulejo&amp;utm_campaign=Badge_Grade
    :alt: Codacy.io grade

.. |coverage| image:: https://codecov.io/gh/legumeinfo/azulejo/branch/master/graph/badge.svg
    :target: https://codecov.io/gh/legumeinfo/azulejo
    :alt: Codecov.io test coverage

.. |issues| image:: https://img.shields.io/github/issues/LegumeFederation/lorax.svg
    :target:  https://github.com/legumeinfo/azulejo/issues
    :alt: Issues reported

.. |requires| image:: https://requires.io/github/legumeinfo/azulejo/requirements.svg?branch=master
     :target: https://requires.io/github/legumeinfo/azulejo/requirements/?branch=master
     :alt: Requirements Status

.. |dlrate| image:: https://img.shields.io/pypi/dm/azulejo
    :target: https://pypistats.org/packages/azulejo
    :alt: Download stats

.. |downloads| image:: https://pepy.tech/badge/azulejo
    :target: https://pepy.tech/project/azulejo
    :alt: Download stats
