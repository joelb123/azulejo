#!/bin/bash
#
# Configuration and run script for dagchainer 
#
set -e # stop on errors
version="1.0"
script_name="$(basename "${BASH_SOURCE}")"
script_dir=''
pushd "$(dirname "$(readlink -f "$BASH_SOURCE")")" >/dev/null && {
  script_dir="$PWD"
  popd >/dev/null
}
scriptstart=$(date +%s)
pkg="${script_name%.sh}"
PKG="$(echo ${pkg} | tr /a-z/ /A-Z/)"
PKG_DIR="${PKG}_DIR"
PKG_WORK_DIR="${PKG}_WORK_DIR"
if [ -z "${!PKG_DIR}" ]; then
  root_dir=~/.local/share/${pkg}
else
  root_dir="${!PKG_DIR}"
fi
if [ -z "${!PKG_WORK_DIR}" ]; then
  work_dir="/tmp/${pkg}-work"
else
  work_dir="${!PKG_WORK_DIR}"
fi
etc_dir="${root_dir}/etc"
blast_db_dir="${work_dir}/blast_db"
blast_out_dir="${work_dir}/blast_out"
dag_dir="${work_dir}/dag"
error_exit() {
  echo >&2 "ERROR--unexpected exit from ${BASH_SOURCE} script at line:"
  echo >&2 "   $BASH_COMMAND"
}
trap error_exit EXIT
TOP_DOC="""Compute synteny among sets of GFF/FASTA files

Usage:
        ${pkg} COMMAND [COMMAND_OPTIONS]

Commands (in order they are usually run):
            version - Get installed package version
               init - Initialize parameters required for run
             config - View/set run parameters
                run - Run one or all analysis steps, '-h' to list
       clear_config - Clear all config variables
              clean - Delete work directory

Variables (accessed by \"config\" command):
      blast_threads - threads to use in searches [default: 4]
    dagchainer_args - Argument for DAGchainer command
             dbtype - Database type, either 'nucl' or 'prot'
              e_val - Maximum BLAST score permitted in matches
          fasta_ext - Extension of FASTA files
            gff_ext - Extension of GFF files
            out_dir - output directory [default: './${pkg}_out']
       pct_identity - Minimum percent sequence identity
            version - version of this script at config time

Environmental variables (may be set externally):
         ${PKG}_DIR - Location of the config directory, currently
                       \"${root_dir}\"
     ${PKG}_WORK_DIR - Location of working files, currently
                       \"${work_dir}\"
"""
#
# Helper functions begin here
#
set_value() {
  if [ "$2" == "-d" ]; then
    rm -f "${etc_dir}/${1}"
  else
    echo "$2" >"${etc_dir}/${1}"
  fi
}
get_value() {
  if [ -e ${etc_dir}/${1} ]; then
    cat ${etc_dir}/${1}
  else
    trap - EXIT
    echo >&2 "ERROR--value for $1 variable not found."
    exit 1
  fi
}
howmany() {
  set -f
  set -- $1
  echo $#
}
histogram() {
  if [ $# -lt 1 ]; then
    echo "Usage: histogram  [field] [precision]"
    return 0
  fi

  field=1
  if [ $# -gt 1 ]; then
    field=$2
  fi

  printf '#N\tSize\n'
  cat $1 | sort -k $field -n -t ',' | cut -d ',' -f $field | uniq -c
}
#
perl_defs() {
  #
  # Perl local::lib settings, where Bioperl::SeqIO is installed
  #
  perlbase=${SYN_BIN}/perl5
  PATH="${perlbase}/bin${PATH:+:${PATH}}"; export PATH
  PERL5LIB="${perlbase}/lib/perl5${PERL5LIB:+:${PERL5LIB}}"; export PERL5LIB
  PERL_LOCAL_LIB_ROOT="${perlbase}${PERL_LOCAL_LIB_ROOT:+:${PERL_LOCAL_LIB_ROOT}}"
  export PERL_LOCAL_LIB_ROOT
  PERL_MB_OPT="--install_base \"${perlbase}\""; export PERL_MB_OPT
  PERL_MM_OPT="INSTALL_BASE=${perlbase}"; export PERL_MM_OPT
}
#
# run functions
#
run_ingest() {
  gff_ext=$(get_value gff_ext)
  fasta_ext=$(get_value fasta_ext)
  gff_files="$(ls *.${gff_ext})"
  fasta_files="$(ls *.${fasta_ext})"
  n_fasta=$(howmany "$fasta_files")
  echo "ingest--combining info from ${n_gff}  ${gff_ext} and ${n_fasta} ${fasta_ext} files"
  for path in $gff_files; do
    base=$(basename $path .${gff_ext})
    base_no_ann=$(echo $base | perl -pe 's/\.ann\d+\.\w+//')
    cat $path | awk -v OFS="\t" '$3=="mRNA" {print $1, $4, $5, $9}' |
      perl -pe 's/ID=([^;]+);.+/$1/' >${work_dir}/${base_no_ann}.bed
  done
  echo "adding positional information to FASTA ids"
  for path in ${work_dir}/*.bed; do
    base=$(basename $path .bed)
    cat $path | awk '{print $4 "\t" $1 "__" $4 "__" $2 "__" $3 }' \
      >${work_dir}/${base}.hsh
    hash_into_fasta_id.pl\
      -fasta ${base}.${fasta_ext}\
      -hash ${work_dir}/${base}.hsh \
      -suff_regex \
      >${work_dir}/${base}.${fasta_ext}
  done
  echo
}
#
run_blast_dbs() {
  echo "blast_dbs--creating BLAST databases"
  start_time=$(date +%s)
  fasta_ext=$(get_value fasta_ext)
  for path in ${work_dir}/*.${fasta_ext}; do
    base=$(basename $path .${fasta_ext})
    makeblastdb -in $path -dbtype $(get_value dbtype) \
      -hash_index -parse_seqids -title \
      $base -out ${blast_db_dir}/$base 1>/dev/null &
  done
  wait
  end_time=$(date +%s)
  set_value db_time_s $((end_time-start_time))
}
#
run_blastall() {
  e_val=$(get_value e_val)
  pct_id=$(get_value pct_identity)
  n_threads=$(get_value blast_threads)
  echo "blastall--doing half-diagonal BLAST at ${e_val} E-value and ${pct_id}% identity using ${n_threads} threads"
  start_time=$(date +%s)
  fasta_ext=$(get_value fasta_ext)
  for qry_path in ${work_dir}/*.${fasta_ext}; do
    qry_base=$(basename $qry_path .${fasta_ext})
    for sbj_path in ${work_dir}/*.${fasta_ext}; do
      sbj_base=$(basename $sbj_path .${fasta_ext})
      if [[ "$qry_base" > "$sbj_base" ]]; then
         blastn -query $qry_path -db ${blast_db_dir}/$sbj_base \
          -perc_identity $pct_id \
          -evalue $e_val \
          -outfmt 6 \
          -num_threads $n_threads \
          -out ${blast_out_dir}/${qry_base}.x.${sbj_base}.bln #&
      fi
    done
  done
  end_time=$(date +%s)
  set_value blast_time_s $((end_time-start_time))
}
#
run_filter(){
  pct_id=$(get_value pct_identity)
  echo "filtering top hits at ${pct_id}% sequence identity"
  for blast_path in ${blast_out_dir}/*; do
    outfilebase=$(basename $blast_path .bln)
    cat $blast_path | \
	    awk -v OFS="\t" '$3>='"${pct_id}"' {print $1, $2, $11}' | \
	    top_blast_hit.awk | \
	    perl -pe 's/__/\t/g' > ${dag_dir}/${outfilebase}_matches.tsv
  done
}
#
run_dagchainer() {
  dag_args=$(get_value dagchainer_args)
  echo "dagchainer--running DAGchainer using args ${dag_args} "
  start_time=$(date +%s)
  for match_path in ${dag_dir}/*_matches.tsv; do
    run_DAG_chainer.pl $dag_args \
      -i $match_path 2>/dev/null 1>/dev/null &
  done
  wait
  end_time=$(date +%s)
  set_value dag_time_s $((end_time-start_time))
  echo "generating single-linkage synteny anchors"
  printf "matches\tscore\trev\tid1\tid2\n" >${work_dir}/synteny_blocks.tsv
  for path in ${dag_dir}/*.aligncoords; do
    file1=${path%%.x.*}
    file2=$(basename ${path##*.x.} _matches.tsv.aligncoords)
    cat $path | awk '$1!~/^#/ {print $2 "\t" $6}' \
      >>${work_dir}/homology_pairs.tsv
    cat $path | grep \#\# | grep -v reverse |
      awk '{print substr($14,0,length($14)-2) "\t" $10 "\t" 1 "\t" $3 "\t" $5}' \
        >>${work_dir}/synteny_blocks.tsv
    cat $path | grep \#\# | grep reverse |
      awk '{print substr($15,0,length($15)-2) "\t" $11 "\t" 1 "\t" $3 "\t" $5}' \
        >>${work_dir}/synteny_blocks.tsv
  done
  blinkPerl_v1.1.pl -in ${work_dir}/homology_pairs.tsv \
    -sum ${work_dir}/cluster_sizes.txt \
    -out ${work_dir}/clusters.txt
  echo
}
#
run_summarize() {
  echo "summarize--formatting synteny data and computing stats"
  out_dir="$(get_value out_dir)"
  if [ ! -d "$out_dir" ]; then
      echo "creating output directory \"${out_dir}/\""
      mkdir -p $out_dir
  fi
  rename_reorder_adjacency.py ${work_dir}/clusters.txt ${out_dir}/synteny_anchors.tsv
  n_anchors=$(cat ${out_dir}/synteny_anchors.tsv | wc -l)
  cp ${work_dir}/synteny_blocks.tsv ${out_dir}/synteny_blocks.tsv
  n_blocks=$(cat ${out_dir}/synteny_blocks.tsv | wc -l)
  echo "${n_anchors} anchors produced ${n_blocks} synteny blocks"
  echo "Computing homology clusters"
  pairs_to_adjacency.py ${work_dir}/homology_pairs.tsv ${out_dir}/homology_clusters.tsv
  join_homology_synteny.py ${work_dir}/synteny_anchors.tsv ${work_dir}/homology_clusters.tsv
  printf "stat\tvalue\n" >${out_dir}/stats.tsv
  head -2 ${work_dir}/cluster_sizes.txt |
    sed -e 's/: /\t/' |
    sed -e 's/\.//' |
    sed -e 's/Total number of //' \
      >>${out_dir}/stats.tsv
  seqids=$(grep seqids ${work_dir}/cluster_sizes.txt | awk '{print $8}')
  scriptend=$(date +%s)
  printf "seqids_in_clusters\t$seqids\n" >>${out_dir}/stats.tsv
  printf "db_time_s\t$(get_value db_time_s)\n" >>${out_dir}/stats.tsv
  printf "blast_time_s\t$(get_value blast_time_s)\n" >>${out_dir}/stats.tsv
  printf "dag_time_s\t$(get_value dag_time_s)\n" >>${out_dir}/stats.tsv

  echo
  cat ${out_dir}/stats.tsv
}
#
# top-level command functions
#
config() {
  CONFIG_DOC="""Sets/displays key/value pairs for the $pkg build system.

Usage:
   $scriptname set [-h] [KEY] [VALUE | -d]

Options:
   -h   Prints this help message and exits
   -d   Deletes the setting of KEY


Arguments:
   If KEY is absent, all values will be displayedWhat
   If KEY is present but VALUE is absent, the value will be displayed
   If KEY and VALUE are present, the value will be set
"""
  if [ "$#" -eq 0 ]; then
    echo >&2 "$CONFIG_DOC"
    param="all"
  elif [ "$1" == "-h" ]; then
    echo >&2 "$CONFIG_DOC"
    param=all
  else
    param="$1"
  fi
  if [ "$param" == "all" ]; then
      trap - EXIT
      for key in $(ls ${etc_dir}); do
        value="$(get_value ${key})"
        printf '%-20s\t%s\n' ${key} ${value} >&1
      done
      exit 0
  fi
  if [ "$#" -eq 1 ]; then
    if [ -e ${etc_dir}/${param} ]; then
      echo "$(get_value $param)"
    else
      trap - EXIT
      echo >&2 "ERROR--\"${1}\" has not been set"
      exit 1
    fi
  elif [ "$#" -eq 2 ]; then # set
    set_value $param $2
  else
    trap - EXIT
    echo >&2 "$CONFIG_DOC"
    echo >&2 "ERROR--too many arguments (${#})."
    exit 1
  fi
}
#
clear_config() {
  echo "clearing configuration directory"
  rm -f ${etc_dir}/*
}
#
init() {
  echo "setting run configuration parameters"
  echo
  set_value e_val "1e-10"
  set_value blast_threads 4
  set_value pct_identity 95
  set_value dagchainer_args ""
  set_value fasta_ext "fna"
  set_value gff_ext "gff3"
  set_value dbtype "nucl"
  set_value out_dir "${pkg}_out"
  config all
}
#
run() {
  RUN_DOC="""Run an analysis step

Usage:
   $scriptname run [STEP]

Steps:
   If STEP is not set, the following steps will be run in order,
   otherwise the step is run by itself:
              ingest - get info from matching GFF and FASTA files
           blast_dbs - create BLAST databases
            blastall - do all-against-all blast
              filter - reduce to top reciprocal hits above pct_identity
          dagchainer - compute Directed Acyclic Graphs
           summarize - compute synteny stats
"""
  commandlist="ingest blast_dbs blastall filter dagchainer summarize"
  if [ "$#" -eq 0 ]; then # run the whole list
    for package in $commandlist; do
      run_$package
      echo
    done
  else
    command="$1"
    shift 1
    case $commandlist in
    *"$command"*)
      run_$command $@
      ;;
    $commandlist)
      trap - EXIT
      echo >&2 "$RUN_DOC"
      if [ "$command" == "-h" ]; then
        exit 0
      fi
      echo >&2 "ERROR--unrecognized run step \"$1\""
      exit 1
      ;;
    esac
  fi
}
#
version() {
  echo $version
}
#
clean() {
  echo "cleaning work directory and tmpOut files"
  rm -rf $work_dir 
  rm -f .*.tmpOut
}
#
# Command-line interpreter.
#
if [ "$#" -eq 0 ]; then
  trap - EXIT
  echo >&2 "$TOP_DOC"
  exit 1
fi
# Create directories if needed
dirlist="root_dir work_dir etc_dir blast_db_dir blast_out_dir dag_dir"
for dirvar in $dirlist; do
    dirname="${!dirvar}"
    if [ ! -d "$dirname" ]; then
      echo "creating directory \"${dirname}\" as $dirvar"
      mkdir -p $dirname
    fi
done
#
command="$1"
shift 1
case $command in
"config")
  config $@
  ;;
"clean")
  clean $@
  ;;
"clear_config")
  clear_config $@
  ;;
"init")
  init $@
  ;;
"run")
  run $@
  ;;
"version")
  version $@
  ;;
*)
  trap - EXIT
  echo >&2 "ERROR--command \"$command\" not recognized."
  exit 1
  ;;
esac
trap - EXIT
exit 0
