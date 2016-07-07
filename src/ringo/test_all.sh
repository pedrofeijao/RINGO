#!/usr/bin/env bash
# Usage:
#   test_all.sh
#
# Depends on:
#  Optionally:
#  - MGRA
#  - DeClone
#
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
_indel_perc=0.2
_scale=1.5
_disturb=1
_indel_length=5
_n_chr=1
_ringo_weights="custom_weight.txt"
_kt=0.4

###############################################################################
# Available software:
###############################################################################
_declone=$(awk -F "=" '/declone_path/ {print $2}' ringo.cfg)
_mgra=$(awk -F "=" '/mgra_path/ {print $2}' ringo.cfg)
###############################################################################
# Options from config:
###############################################################################
_declone_weights=$(awk -F "=" '/declone_internal_weight/ {print $2}' ringo.cfg)

###############################################################################
# Program Functions
###############################################################################
_run_all() {

  _simulation

  _ringo_weights

  _declone_weights

  _run_ringo

  _run_ringo_declone

  _run_mgra

}

_simulation() {
  printf "Simulating  a dataset...\n"
  ./simulation.py -n $_n_genes -c $_n_chr -o $_out_folder -s $_n_genomes -i $_indel_perc -sc $_scale -d $_disturb --indel_length $_indel_length
  printf "Done.\n\n"
}

_ringo_weights() {
  printf "Generating RINGO weights ...\n"
  ./ancestral_weights.py -i $_out_folder/leaf_genomes.txt -t $_out_folder/evolved_tree.nwk -o $_out_folder/$_ringo_weights
  printf "Done.\n\n"
}

_declone_weights() {
  printf "Trying to run DeClone...\n"
  if command -v $_declone >/dev/null 2>&1; then
    printf "DeClone found, generating DeClone weights ...\n"
    ./weighting_DeClone.py -t $_out_folder/evolved_tree.nwk -m $_out_folder/leaf_genomes.txt -kT $_kt -o $_out_folder
    printf "Done.\n\n"
  else
    printf "DeClone not found, skipping. If you want to run DeClone, please install it and edit 'ringo.cfg' to indicate the executable path.\n\n"
  fi

}

_run_ringo() {
  printf "Running RINGO on default weights ...\n"
  ./ringo.py -i $_out_folder/leaf_genomes.txt -t $_out_folder/evolved_tree.nwk -o $_out_folder -w $_out_folder/$_ringo_weights
  printf "Done.\n\n"
}

_run_ringo_declone() {
  if command -v $_declone >/dev/null 2>&1; then
    printf "Running RINGO on DeClone weights ...\n"
    ./ringo.py -i $_out_folder/leaf_genomes.txt -t $_out_folder/evolved_tree.nwk -o $_out_folder -w $_out_folder/${_declone_weights}${_kt}
    printf "Done.\n\n"
  fi
}

_run_mgra() {
  printf "Trying to run MGRA...\n"
  if command -v $_mgra >/dev/null 2>&1; then
    printf "Running MGRA...\n"
    printf "Done.\n\n"
  else
    printf "MGRA not found, skipping. If you want to run MGRA, please install it and edit 'ringo.cfg' to indicate the executable path.\n\n"
  fi
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
