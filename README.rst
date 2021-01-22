.. epigraph:: azulejo
              noun INFORMAL
              a glazed tile, usually blue, found on the inside of churches and palaces in Spain and Portugal.

azulejo
=======
``azulejo`` azulejo combines homology and synteny information to
tile phylogenetic space.
The inputs to ``azulejo`` are FASTA files of nucleotide-space
sequences of primary-transcript protein genes and their associated GFF files.
Outputs are sets of proxy gene fragments chosen for 
their concordance in multiple sequence alignments, along with
subtrees.

Prerequisites
-------------
Python 3.7 or greater is required. ``azulejo`` is tested under Linux 
using Python 3.8 and 3.9 and under MacOS Big Sur using XCode command-line tools
system Python (currently 3.8). Mac users should see the `instructions
on configuring their systems <macos.rst>`_.  Installation on BSD is not
supported because many of the python dependencies lack BSD wheels.

We recommend you install ``azulejo`` into its own virtual environment due
to the large number of python dependencies.  The easiest way for most 
users to install and maintain up-to-date virtual environments is via the
tool `pipx<https://pipxproject.github.io/pipx/>`_.  If your system does
not have ``pipx`` installed, you can do so via the commands::

        python3 -m pip install --user --upgrade pip
        python3 -m pip install --user --upgrade pipx
        python3 -m pipx ensurepath

Follow any instructions that the last command produces about starting a new
shell if necessary.  

If you choose to have ``azulejo`` compile and install its binary dependencies,
you will need compilers, ``make``, and ``cmake``  and standard headers
for ``zlib`` and ``bz2``.  All linux systems configured for development will have
these available.  We test compilation under gcc version 10.2 on linux and
clang 12.0.0 on MacOS.  We use program-guided optimization for one of the
binary dependencies, and we believe that gcc 10 does a much better job
of optimization than gcc 9, so it may benefit you to upgrade your compiler
if needed.

Installation for Users
----------------------
Once the prerequisite has been met, you may then install ``azulejo`` 
in its own virtual environment by issuing the command::

        pipx install azulejo

``azulejo`` contains some long commands and many options.  To enable command-line
completion for ``azulejo`` commands, execute the following command if you are using
``bash`` as your shell: ::

    eval "$(_AZULEJO_COMPLETE=source_bash azulejo)"

Then you should run ``azulejo install`` and check the versions of all binary
dependencies that may installed system-wide.


Environmental Variables
-----------------------
``azulejo`` recognizes the following environmental variables:

* AZULEJO_INSTALL_DIR
This is a writable directory for installation of binary dependencies.  Binaries
will go into the ``bin`` directory.  The default is the virtual environment
directory.

* BUILD_DEV
This is the directory used for building binary dependencies.  Default is the
first memory device found for linux (e.g., `/run/shm`) or `/tmp` for MacOS.
Set this if compilation fails because it runs out of memory.

* SCRATCH_DEV
This is the directory used for temporary merging of lists.  The default is
`/tmp`, but you may set it to a fast memory based device if you have enough
memory.

* MAKEOPTS
These are the arguments to the ``make`` and ``make install`` commands when
building dependencies.  It's good to set this to the number of processors
on your system via the command ``export MAKEOPTS="-j $(nproc)"`` to speed
up installation.  The only time this variable is used is during
``azulejo install``.

* SPINNER_UPDATE_PERIOD
This is the number of seconds between updates of the spinner.  This
defaults to 1, but it is advisable to set it higher for automated testing
so as not to exceed logfile character limits.

* LOG_TO_PRINT
If set, the logger will be a simple print without using the more
complex functions of ``loguru`` such as colors and logging to files.
This is sometimes useful in automated testing.



Installation of Binary Dependencies
-----------------------------------
``azulejo`` requires `MMseqs <https://github.com/soedinglab/MMseqs2>`_ 
for homology clustering and `MUSCLE <https://www.drive5.com/muscle/downloads.htm>`_
for sequence alignment and initial tree-building.
``azulejo`` installs binaries into the virtualenv by default, so
any systemwide installations of these packages will not get clobbered by the install.
In particular, ``muscle`` is PGO-optimized, which gives nearly a factor of 2 higher
performance than prebuilt binaries.  We recommand you set ``MAKEOPTS`` as explained
above, then issue the command ``azulejo install all`` to ensure you get correct versions
optimized for your hardware.


There are three optional dependencies that can be installed via ``azulejo install`` 
that are of interest only to a small subset of users who wish to compare against
other homology clustering and synteny methods.  
`usearch <https://www.drive5.com/usearch/download.html>`_ 
is a licensed homology clustering program that is free for individual, non-commercial
use that can be downloaded and installed by the ``azulejo install usearch``
command after accepting the license terms.  ``azulejo install dagchainer-tool`` gets you
a somewhat crude Bash script that uses BLAST homology clustering followed by 
synteny calculation via `DAGchainer <https://dagchainer.sourceforge.net>`_.  
``dagchainer-tool`` will need the dependency of ``perl`` with ``bioperl`` installed.
``dagchainer_tool`` increases the sequence ID length as part of its processing, so
if any of your sequence IDS are longer than about 30 characters, they will violate BLAST's
hard limit of 50 characters in sequence ID fields.  In that case you will need
to install a patched version of BLAST using the command ``azulejo install blast-longids``.

Installation For Developers
---------------------------
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
