#!/bin/bash

# Prune haplotype files used for calculating power

pop=AFR
#pruning=pruneSepRaresim
nsim=20000
pcase=160
pexp=100
pconf=80

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
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pcase}.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pexp}.txt \
    -l chr19.block37.${pop}.sim${rep}.copy.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pexp}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.gz
    #-L ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pexp}syn.legend \
    #-H ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.gz

# Convert -H hap file to sm
python3 ${WD}/raresim/convert.py \
    -i ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.gz \
    -o ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.sm
    #-i ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.gz \
    #-o ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.sm 

# Prune fun variants down to pexp %
python3 ${WD}/raresim/sim.py \
    -m ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.sm \
    --f_only ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pexp}_6bins.txt \
    -l ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pexp}syn.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pexp}fun.${pexp}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.gz
    #-m ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pexp}syn.haps.sm \
    #-l ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pexp}syn.legend \
    #-L ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pexp}fun.${pexp}syn.legend \
    #-H ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.gz

# Convert -H hap file to sm
python3 ${WD}/raresim/convert.py \
    -i ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.gz \
    -o ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.sm 
    #-i ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.gz \
    #-o ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.sm 

# Prune fun and syn variants down again to pconf %
python3 ${WD}/raresim/sim.py \
    -m ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pconf}_6bins.txt \
    -l ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pexp}fun.${pexp}syn.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz
    #-m ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pexp}fun.${pexp}syn.haps.sm \
    #-l ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pexp}fun.${pexp}syn.legend \
    #-L ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    #-H ${WD}/mastersProject/20K_${pop}/${pruning}/${pcase}v${pexp}v${pconf}/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz

done
