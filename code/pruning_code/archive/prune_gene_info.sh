#!/bin/bash

# Code to prune hap files to 100% fun 100% syn for Block 37 Gene Info purposes

pop=NFE
#pruning=pruneSepRaresim
nsim=20000
#pcase=120
pexp=100
#pconf=80

# define parameters
end=100
start=$(($end-99))

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa

cd /storage/math/projects/RAREsim/Cases/Sim_20k/${pop}/data

### use RAREsim to prune to pfun % functional and psyn % synonymous
for rep in $(eval echo "{$start..$end}")
do

#/storage/singularity/mixtures.sif
# For pruning Separately using only RAREsim-python
# Prune fun variants down to pcase % and syn variants down to pexp %
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop}.sim${rep}.controls.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pexp}.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pexp}.txt \
    -l chr19.block37.${pop}.sim${rep}.copy.legend \
    -L ${WD}/mastersProject/gene_info/20K_${pop}/chr19.block37.${pop}.sim${rep}.${pexp}fun.${pexp}syn.legend \
    -H ${WD}/mastersProject/gene_info/20K_${pop}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.gz

done
