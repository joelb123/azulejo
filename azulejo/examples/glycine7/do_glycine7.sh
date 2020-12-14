#!/bin/bash
set -e
echo "This script reads 7 glycine genomes from the Legume Federation Data Store"
echo "and calculates homology and synteny for them."
echo "installing binary dependencies, if needed"
azulejo install all
echo "finding the files in the data store"
azulejo find-files -po  -pr genus site://legfed/Glycine_max glycine7 glycine7.toml
azulejo find-files site://legfed/Glycine_soja glycine7.glyso glycine7.toml
azulejo find-files -e '*Wm82.gnm2.ann2*' -e '*Lee*' -e '*Zh13*' -pr strain -r minstrain \
	-n "g4" -n "g2a1" -n "g1" site://legfed/Glycine_max glycine7.glyma.Wm82 glycine7.toml
echo "reading the data into local files"
azulejo ingest glycine7.toml
echo "calculating homology info, including alignments and crude trees"
azulejo homology glycine7
echo "doing the synteny calculation"
azulejo synteny glycine7
echo "Done with glycine7"
