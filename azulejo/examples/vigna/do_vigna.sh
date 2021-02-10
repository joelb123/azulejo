#!/bin/bash
set -e
echo "This script reads vigna genomes from the Legume Federation Data Store"
echo "and calculates homology and synteny for them."
azulejo install all
azulejo find-files -po  -pr genus site://legfed/Vigna_angularis vigna vigna.toml
azulejo find-files site://legfed/Vigna_angularis vigna.angularis vigna.toml
azulejo find-files -a '*.protein.faa.gz' site://legfed/Vigna_radiata vigna.radiata vigna.toml
azulejo find-files site://legfed/Vigna_unguiculata vigna.unguiculata vigna.toml
azulejo ingest vigna.toml
azulejo homology vigna
azulejo synteny vigna
echo "Done with vigna"
