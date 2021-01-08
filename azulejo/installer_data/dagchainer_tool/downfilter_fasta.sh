#!/bin/bash
script=$(basename $0)
DOC="""${script}--Downselect FASTA records to match those in another file

USAGE:
   ${script} IDFILE INFILE OUTFILE

   where
     IDFILE is a FASTA file whose IDS are the ones you want
     INFILE is a FASTA file that will be downselected
    OUTFILE is the downselected FASTA file

DEPENDENCIES:
     fastaq (python package)"""
if [ "$#" -lt 3 ]; then
  echo >&2 "$DOC"
  exit 1
fi
set -e
IDFILE=$1
INFILE=$2
OUTFILE=$3
LONG_IDS=$(mktemp)
IDS=$(mktemp)
fastaq get_ids $IDFILE $LONG_IDS
awk '{print $1}' $LONG_IDS > $IDS
fastaq filter --ids_file $IDS $DOWNSELECTFILE $OUTFILE
incount=$(grep \> $INFILE | wc -l)
outcount=$(grep \> $OUTFILE | wc -l)
rm -f $LONG_IDS $IDS
echo "${incount} FASTA records in, ${outcount} out"
exit 0
