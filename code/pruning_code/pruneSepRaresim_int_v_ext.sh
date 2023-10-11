#!/bin/bash

# Test 99% internal v 99% external using internal and external samples from 2 separate 
# 99% pruned hap files, but both stemmin from same 100% pruned hap file

pop=NFE
pruning=pruneSepRaresim
nsim=20000
pcase=100
pconf=99
int_prune=100
ext_prune=99


# define parameters
end=100
start=$(($end-99))

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa

cd /storage/math/projects/RAREsim/Cases/Sim_20k/${pop}/data

### Loop through all simulation replicates
for rep in $(eval echo "{$start..$end}")
do

#/storage/singularity/mixtures.sif
# For pruning Separately using only RAREsim-python
# Prune fun and syn variants down to pcase %
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop}.sim${rep}.controls.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pcase}.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pcase}.txt \
    -l chr19.block37.${pop}.sim${rep}.copy.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz

# Convert -H hap file to sm
python3 ${WD}/raresim/convert.py \
    -i ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz \
    -o ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm 

# Prune fun and syn variants down again to pconf % for INTERNAL samples
python3 ${WD}/raresim/sim.py \
    -m ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pconf}_6bins.txt \
    -l ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/internal_data/chr19.block37.${pop}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/internal_data/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz

# Prune fun and syn variants down again to pconf % for EXTERNAL samples
python3 ${WD}/raresim/sim.py \
    -m ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pconf}_6bins.txt \
    -l ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/external_data/chr19.block37.${pop}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pruning}/int_v_ext/external_data/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz

done
