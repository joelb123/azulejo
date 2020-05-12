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
Python 3.6 or greater is required.
This package is tested under Linux and MacOS using Python 3.7.

Installation for Users
----------------------
Install via pip or (better yet) `pipx <https://pipxproject.github.io/pipx/>`_: ::

     pipx install azulejo

``azulejo`` contains some long commands and many options.  To enable command-line
completion for ``azulejo`` commands, execute the following command if you are using
``bash`` as your shell: ::

    eval "$(_AZULEJO_COMPLETE=source_bash azulejo)"

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
The currently implemented commands are:

========================= ==================================================
  add-singletons          Add singleton clusters to cluster file.
  adjacency-to-graph      Turn adjacency data into GML graph file.
  analyze-clusters        Statistics of clustering as function of identity.
  annotate-homology       Marshal homology and sequence information.
  cluster-in-steps        Cluster in steps from low to 100% identity.
  clusters-to-histograms  Compute histograms from cluster file.
  combine-clusters        Combine synteny and homology clusters,
  compare-clusters        compare one cluster file with another
  dagchainer-synteny      Read DAGchainer synteny into homology frames.
  length-std-dist         Plot length distribution of singletons in clusters
  outlier-length-dist     Plot length distribution of outliers in clusters.
  prepare-protein-files   Sanitize and combine protein FASTA files.
  proxy-genes             Calculate a set of proxy genes from synteny files.
  synteny-anchors         Calculate synteny anchors.
  usearch-cluster         Cluster at a global sequence identity threshold.
========================= ==================================================

Each command has its ``COMMANDOPTIONS``, which may be listed with: ::

    azulejo COMMAND --help

Project Status
--------------
+-------------------+------------+------------+
| Latest Release    | |pypi|     | |azulejo|  |
+-------------------+------------+            +
| GitHub            | |repo|     |            |
+-------------------+------------+            +
| License           | |license|  |            |
+-------------------+------------+            +
| Travis Build      | |travis|   |            |
+-------------------+------------+            +
| Coverage          | |coverage| |            |
+-------------------+------------+            +
| Code Grade        | |codacy|   |            |
+-------------------+------------+            +
| Dependencies      | |depend|   |            |
+-------------------+------------+            +
| Issues            | |issues|   |            |
+-------------------+------------+------------+


.. |azulejo| image:: docs/azulejo.jpg
     :target: https://en.wikipedia.org/wiki/Azulejo
     :alt: azulejo Definition

.. |pypi| image:: https://img.shields.io/pypi/v/azulejo.svg
    :target: https://pypi.python.org/pypi/azulejo
    :alt: Python package

.. |repo| image:: https://img.shields.io/github/commits-since/legumeinfo/azulejo/0.3.svg
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

.. |depend| image:: https://api.dependabot.com/badges/status?host=github&repo=legumeinfo/azulejo
     :target: https://app.dependabot.com/accounts/legumeinfo/repos/203668510
     :alt: dependabot dependencies
