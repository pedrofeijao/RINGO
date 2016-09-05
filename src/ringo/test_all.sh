#!/usr/bin/env bash
# Usage:
#   test_all.sh
#
# Depends on:
#  Optionally:
#  - MGRA - https://github.com/ablab/mgra
#  - DeClone - https://github.com/yannponty/DeClone
#  - PhySca - https://github.com/nluhmann/PhySca

# Bash Boilerplate: https://github.com/alphabetum/bash-boilerplate
#
# Copyright (c) 2015 William Melody • hi@williammelody.com


###############################################################################
# Environment
###############################################################################

# $_ME
#
# Set to the program's basename.
_ME=$(basename "${0}")

###############################################################################
# Help
###############################################################################

_print_help() {
  cat <<HEREDOC
Usage:
  ${_ME} [<arguments>]
  ${_ME} -h | --help

Options:
  -h --help  Show this screen.
HEREDOC
}
###############################################################################
# Default Parameters
###############################################################################
_out_folder="test_all_output"
_n_genomes=6
_n_genes=100
_n_chr=1
_ins_perc=0.1
_del_perc=0.1
_scale=1.5
_disturb=1
_indel_length=4
_n_chr=1
_ringo_weights="custom_weight.txt"
_kt=0.4

###############################################################################
# Available software:
###############################################################################
_declone=$(awk -F "=" '/declone_path/ {print $2}' ringo.cfg)
_mgra=$(awk -F "=" '/mgra_path/ {print $2}' ringo.cfg)
_physca=$(awk -F "=" '/physca_path/ {print $2}' ringo.cfg)

###############################################################################
# Options from config:
###############################################################################
_declone_weights=$(awk -F "=" '/declone_internal_weight=/ {print $2}' ringo.cfg)
_declone_extant_weights=$(awk -F "=" '/declone_extant_weight=/ {print $2}' ringo.cfg)
_sim_tree=$(awk -F "=" '/sim_tree=/ {print $2}' ringo.cfg)
_extant_genomes=$(awk -F "=" '/sim_extant_genomes=/ {print $2}' ringo.cfg)
_ancestral_genomes=$(awk -F "=" '/sim_ancestral_genomes=/ {print $2}' ringo.cfg)
_mgra_config=$(awk -F "=" '/sim_mgra_config=/ {print $2}' ringo.cfg)
_mgra_output=$(awk -F "=" '/mgra_output_folder=/ {print $2}' ringo.cfg)
_scj_genomes=$(awk -F "=" '/scj_genomes=/ {print $2}' ringo.cfg)
_physca_output=$(awk -F "=" '/physca_output_folder=/ {print $2}' ringo.cfg)

###############################################################################
# Program Functions
###############################################################################
_run_all() {

  # get input:
  _n_genomes=${1:-$_n_genomes}
  _n_genes=${2:-$_n_genes}

  # output locations:
  _ringo_custom_out="ringo_custom"
  _ringo_declone_out="ringo_declone"

  # store algorithms for parse:
  ALGS=()

  # START:
  _simulation

  _ringo_weights

  _declone_weights

  _run_ringo

  _run_ringo_declone

  _run_scj

  _run_mgra

  _run_physca

  _parse_sim

}

_simulation() {
  printf "=== Simulating a dataset...\n"
  ./simulation.py -n $_n_genes -c $_n_chr -o $_out_folder -s $_n_genomes -ip $_ins_perc -dp $_del_perc -sc $_scale -d $_disturb --indel_length $_indel_length
  printf "Done.\n\n"
}

_ringo_weights() {
  printf "=== Generating RINGO weights ...\n"
  ./ancestral_weights.py -i $_out_folder/$_extant_genomes -t $_out_folder/$_sim_tree -o $_out_folder/$_ringo_weights
  printf "Done.\n\n"
}

_declone_weights() {
  printf "=== Trying to run DeClone...\n"
  if command -v $_declone >/dev/null 2>&1; then
    printf "DeClone found, generating DeClone weights ...\n"
    ./weighting_DeClone.py -t $_out_folder/$_sim_tree -m $_out_folder/$_extant_genomes -kT $_kt -o $_out_folder
    printf "Done.\n\n"
  else
    printf "DeClone not found, skipping. If you want to run DeClone, please install it and edit 'ringo.cfg' to indicate the executable path.\n\n"
  fi

}

_run_ringo() {
  printf "=== Running RINGO on default weights ...\n"
  ./ringo.py -i $_out_folder/$_extant_genomes -t $_out_folder/$_sim_tree -o $_out_folder/$_ringo_custom_out -w $_out_folder/$_ringo_weights
  ALGS+=("RINGO_DEFAULT,ringo,$_ringo_custom_out")
  printf "Done.\n\n"
}

_run_ringo_declone() {
  if command -v $_declone >/dev/null 2>&1; then
    printf "=== Running RINGO on DeClone weights ...\n"
    ./ringo.py -i $_out_folder/$_extant_genomes -t $_out_folder/$_sim_tree -o $_out_folder/$_ringo_declone_out -w $_out_folder/${_declone_weights}${_kt}
    ALGS+=("RINGO_DECLONE,ringo,$_ringo_declone_out")
    printf "Done.\n\n"
  fi
}

_run_scj() {
  printf "=== Running SCJ Small Phylogeny...\n"
  ./run_scj.py -i $_out_folder/$_extant_genomes -t $_out_folder/$_sim_tree -o $_out_folder
  ALGS+=("SCJ,scj,$_scj_genomes")
  printf "Done.\n\n"

}
_run_mgra() {
  printf "=== Trying to run MGRA...\n"
  if command -v $_mgra >/dev/null 2>&1; then
    printf "MGRA found, running...\n"

    mkdir -p $_out_folder/$_mgra_output
    command="$_mgra -g $_out_folder/$_extant_genomes -o $_out_folder/$_mgra_output -c $_out_folder/$_mgra_config"
    printf "$command\n"
    $command > $_out_folder/$_mgra_output/mgra2.out
    ALGS+=("MGRA,mgra,$_mgra_output")
    printf "Done.\n\n"
  else
    printf "MGRA not found, skipping. If you want to run MGRA, please install it and edit 'ringo.cfg' to indicate the executable path.\n\n"
  fi
}

_run_physca() {
  printf "=== Trying to run PhySCA...\n"
    if [ -e $_physca ]; then
      if command -v $_declone >/dev/null 2>&1; then
        printf "PhySca found, running...\n"
        mkdir -p $_out_folder/$_physca_output
        command="python $_physca -tree $_out_folder/$_sim_tree -alpha 0.5 -internal $_out_folder/${_declone_weights}${_kt} -extant $_out_folder/${_declone_extant_weights}${_kt} -out $_out_folder/$_physca_output"
        printf "$command\n"
        $command > $_out_folder/$_physca_output/physca.out
        ALGS+=("PhySca,physca,$_physca_output")
        printf "Done.\n\n"
      else
        printf "PhySca found, but DeClone not, it is needed for PhySca. Please install it and edit 'ringo.cfg' to indicate the executable path.\n\n"
      fi
    else
      printf "PhySca not found, skipping. If you want to run PhySca, please install it and edit 'ringo.cfg' to indicate the executable path.\n\n"
    fi

}

_parse_sim() {
  ./parse_sims.py -a ${ALGS[*]} -f test_all_output -t 4
}
###############################################################################
# Main
###############################################################################

# Description:
#   Entry point for the program, handling basic option parsing and dispatching.
_main() {
  # Avoid complex option parsing when only one program option is expected.
  if [[ "${1:-}" =~ ^-h|--help$  ]]
  then
    _print_help
  else
    _run_all "${@}"
  fi
}

# Call `_main` after everything has been defined.
_main "${@:-}"
