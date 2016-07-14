#!/usr/bin/env bash

F=$1 #$1 - folder
S=$2 #$2 - scale
I=$3 #$3 - indel rate
K=0.4 # $4 # declone kT
FILTER=0.05 # $5 # weight filter

F="${F}/sc${S}_indel${I}"

mkdir -p $F

printf "Simulating dataset...\n"
./simulation.py -n 1000 -c 1 -o $F -s 12 -i $I -sc $S -d 1 --indel_length 5

printf "DeClone weights\n"
./weighting_DeClone.py -t $F/sim_tree_no_lengths.nwk -m $F/extant_genomes.txt -kT $K -o $F


printf "Running ringo on simulated branch lenghts...\n"
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree.nwk -o $F/ringo_lengths 2>/dev/null

printf "Running ringo without branch lenghts...\n"
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree_no_lengths.nwk -o $F/ringo_no_lengths 2>/dev/null

printf "Running ringo on branch lenghts estimated with min evol LP...\n"
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree_no_lengths.nwk -o $F/ringo_est_lengths_lp -bl lp 2>/dev/null

printf "Running ringo on branch lenghts estimated with Fitch-Margoliash LSQ...\n"
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree_no_lengths.nwk -o $F/ringo_est_lengths_ls -bl least_squares 2>/dev/null

printf "Running ringo on DeClone weights (branch lengths are not used)...\n"
./ringo.py -i $F/extant_genomes.txt -t $F/sim_tree_no_lengths.nwk -o $F/declone_kt_${K}_f_${FILTER} -w $F/weighted_internal_adjacencies_$K -f $FILTER 2>/dev/null

# ./parse_sims.py -a RINGO,ringo,ringo_lengths RINGO_NO_L,ringo,ringo_no_lengths RINGO_EST,ringo,ringo_est_lengths R_DEC_${K}_f${FILTER},ringo,declone_kt_${K}_f_${FILTER} -f $F  -c
./parse_sims.py -a RINGO,ringo,ringo_lengths RINGO_NO_L,ringo,ringo_no_lengths RINGO_LP,ringo,ringo_est_lengths_lp RINGO_LSQ,ringo,ringo_est_lengths_ls R_DEC_${K}_f${FILTER},ringo,declone_kt_${K}_f_${FILTER} -f $F -t
