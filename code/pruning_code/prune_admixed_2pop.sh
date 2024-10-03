#!/bin/bash

# Code mostly from /storage/math/projects/compinfo/simulations/code/admixed/create_raresim_haps_2pop_v2.sh

# Simulate an admixed African American pop composed of AFR and NFE haplotypes
# Then prune haplotypes accordingly
# Note: Currently need to add -z flag at pruning steps
# extract.py can now handle .gz files (but it can only do .gz files)

pop1=AFR
pop2=NFE
admx_pop1=80
admx_pop2=20
Npop1=13200
Npop2=4800
Nsim_admx=$((${Npop1}+${Npop2}))
scen=s1
sub_scen=default
pcase=160
pconf=80

# define parameters
end=100
start=$(($end-99))

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa

cd ${WD}/admixed/${admx_pop1}${pop1}_${admx_pop2}${pop2}/${scen}/${sub_scen}/pruned_haps

### use RAREsim to prune to pcase % functional and 100% synonymous
for rep in $(eval echo "{$start..$end}")
do 

echo "Simulation Replicate: ${rep}"

### extract the AFR haplotypes
python3 ${WD}/raresim_fix/raresim/extract.py \
    -i ${WD}/Sim_100k/${pop1}/data/chr19.block37.${pop1}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop1}.reduced.sim${rep}.controls.haps.gz \
    -n $((2*$Npop1)) \
    --seed $rep

### extract the NFE haplotypes
python3 ${WD}/raresim_fix/raresim/extract.py \
    -i ${WD}/Sim_100k/${pop2}/data/chr19.block37.${pop2}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop2}.reduced.sim${rep}.controls.haps.gz \
    -n $((2*$Npop2)) \
    --seed $rep

### combine the extracted AFR and NFE haplotype files 
paste <(zcat chr19.block37.${pop1}.reduced.sim${rep}.controls.haps-sample.gz) <(zcat chr19.block37.${pop2}.reduced.sim${rep}.controls.haps-sample.gz) -d " " | gzip > chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz

### convert the combined haplotype file into a sparse matrix
python3 ${WD}/raresim_fix/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.sm

echo "Pruning down to ${pcase}% functional and 100% synonymous variants..."

### prune the haplotypes down to pcase % functional and 100% synonymous variants
# input legend file only needs to change if admixture porportion changes
python3 ${WD}/raresim_fix/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.sm \
    --functional_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_fun_${pcase}_6bins.txt \
    --synonymous_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_syn_100_6bins.txt \
    --stop_threshold 10 \
    -l ${WD}/admixed/${admx_pop1}${pop1}_${admx_pop2}${pop2}/subset_master/chr19.block37.${pop1}_${pop2}.sim${rep}.copy.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.${pcase}fun.100syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz \
    -z

### convert the pruned haplotype file into a sparse matrix
python3 ${WD}/raresim_fix/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.sm

# echo "Pruning down to 100% functional variants..."

# ### prune the haplotypes down to 100% functional
python3 ${WD}/raresim_fix/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.sm \
    --f_only ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_fun_100_6bins.txt \
    --stop_threshold 10 \
    -l chr19.block37.${pop1}_${pop2}.sim${rep}.${pcase}fun.100syn.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.100fun.100syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.gz \
    -z

### convert the pruned haplotype file into a sparse matrix
python3 ${WD}/raresim_fix/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.sm

echo "Pruning down to ${pconf}% functional and ${pconf}% synonymous variants..."

### prune the haplotypes down to pconf % functional and pconf % synonymous variants
python3 ${WD}/raresim_fix/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.sm \
    --functional_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_syn_${pconf}_6bins.txt \
    --stop_threshold 10 \
    -l chr19.block37.${pop1}_${pop2}.sim${rep}.100fun.100syn.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz \
    -z

done

