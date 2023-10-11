#!/bin/bash

# Convert Megan's .hap.gz files to .haps.sm files using RAREsim v2.1.1
# Then prune down to 100% fun and 100% syn
# Then convert to a .sm file
# Then prune down to 80% fun and 80% syn

pop=NFE
pruning=pruneSepRaresim
nsim=20000
pcase=100
pconf=80
int_prune=100
ext_prune=80


# define parameters
end=100
start=$(($end-99))

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa
FD=/home/math/siglersa/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/attempt4_convert_mn_haps_sm
cd /storage/math/projects/RAREsim/Cases/Sim_20k/${pop}/data

### Loop through all replicates
for rep in $(eval echo "{$start..$end}")
do

# Convert Megan's hap file to sm
python3 ${WD}/raresim/convert.py \
    -i chr19.block37.${pop}.sim${rep}.controls.haps.gz \
    -o ${FD}/chr19.block37.${pop}.sim${rep}.controls.haps.sm 

# Prune fun and syn variants down to pcase %
python3 ${WD}/raresim/sim.py \
    -m ${FD}/chr19.block37.${pop}.sim${rep}.controls.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pcase}.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pcase}.txt \
    -l chr19.block37.${pop}.sim${rep}.copy.legend \
    -L ${FD}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -H ${FD}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz 

# Convert -H hap file to sm
python3 ${WD}/raresim/convert.py \
    -i ${FD}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz \
    -o ${FD}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm 

# Prune fun and syn variants down again to pconf %
python3 ${WD}/raresim/sim.py \
    -m ${FD}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pconf}_6bins.txt \
    -l ${FD}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -L ${FD}/chr19.block37.${pop}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H ${FD}/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz

done
