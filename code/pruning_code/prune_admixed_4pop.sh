#!/bin/bash

# Simulate an admixed population composed of 4 different continental ancestry haplotypes
# Then prune haplotypes accordingly
# Note: Currently need to add -z flag at pruning steps
# extract.py can now handle .gz files (but it can only do .gz files)

pop_admx=LTX
pop1=IAM
pop2=NFE
pop3=EAS
pop4=AFR
Npop1=8580
Npop2=8160
Npop3=2700
Npop4=2560
Nsim_admx=$((${Npop1}+${Npop2}+${Npop3}+${Npop4}))
admx_folder=47IAM_44NFE_5EAS_4AFR
scen=s1
sub_scen=default
pcase=160
pconf=80

# define parameters
# end=100
# start=$(($end-99))
rep=1

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa

cd ${WD}/admixed/${admx_folder}/${scen}/${sub_scen}/pruned_haps

### use RAREsim to prune to pcase % functional and 100% synonymous
# for rep in $(eval echo "{$start..$end}")
# do 

echo "Simulation Replicate: ${rep}"

### extract the IAM haplotypes
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

### extract the EAS haplotypes
python3 ${WD}/raresim_fix/raresim/extract.py \
    -i ${WD}/Sim_100k/${pop3}/data/chr19.block37.${pop3}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop3}.reduced.sim${rep}.controls.haps.gz \
    -n $((2*$Npop3)) \
    --seed $rep

### extract the AFR haplotypes
python3 ${WD}/raresim_fix/raresim/extract.py \
    -i ${WD}/Sim_100k/${pop4}/data/chr19.block37.${pop4}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop4}.reduced.sim${rep}.controls.haps.gz \
    -n $((2*$Npop4)) \
    --seed $rep

### combine the extracted haplotype files 
paste <(zcat chr19.block37.${pop1}.reduced.sim${rep}.controls.haps-sample.gz) <(zcat chr19.block37.${pop2}.reduced.sim${rep}.controls.haps-sample.gz) <(zcat chr19.block37.${pop3}.reduced.sim${rep}.controls.haps-sample.gz) <(zcat chr19.block37.${pop4}.reduced.sim${rep}.controls.haps-sample.gz) -d " " | gzip > chr19.block37.${pop_admx}.sim${rep}.controls.haps.gz

### convert the combined haplotype file into a sparse matrix
python3 ${WD}/raresim_fix/raresim/convert.py \
    -i chr19.block37.${pop_admx}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop_admx}.sim${rep}.controls.haps.sm

echo "Pruning down to ${pcase}% functional and 100% synonymous variants..."

### prune the haplotypes down to pcase % functional and 100% synonymous variants
# input legend file only needs to change if admixture porportion changes
python3 ${WD}/raresim_fix/raresim/sim.py \
    -m chr19.block37.${pop_admx}.sim${rep}.controls.haps.sm \
    --functional_bins ${WD}/mac_bin_estimates/${pop_admx}/MAC_bin_estimates_${Nsim_admx}_${pop_admx}_fun_${pcase}_6bins.txt \
    --synonymous_bins ${WD}/mac_bin_estimates/${pop_admx}/MAC_bin_estimates_${Nsim_admx}_${pop_admx}_syn_100_6bins.txt \
    --stop_threshold 10 \
    -l ${WD}/admixed/${admx_folder}/subset_master/chr19.block37.${pop_admx}.sim${rep}.copy.legend \
    -L chr19.block37.${pop_admx}.sim${rep}.${pcase}fun.100syn.legend \
    -H chr19.block37.${pop_admx}.sim${rep}.all.${pcase}fun.100syn.haps.gz \
    -z

### convert the pruned haplotype file into a sparse matrix
python3 ${WD}/raresim_fix/raresim/convert.py \
    -i chr19.block37.${pop_admx}.sim${rep}.all.${pcase}fun.100syn.haps.gz \
    -o chr19.block37.${pop_admx}.sim${rep}.all.${pcase}fun.100syn.haps.sm

echo "Pruning down to 100% functional variants..."

### prune the haplotypes down to 100% functional
python3 ${WD}/raresim_fix/raresim/sim.py \
    -m chr19.block37.${pop_admx}.sim${rep}.all.${pcase}fun.100syn.haps.sm \
    --f_only ${WD}/mac_bin_estimates/${pop_admx}/MAC_bin_estimates_${Nsim_admx}_${pop_admx}_fun_100_6bins.txt \
    --stop_threshold 10 \
    -l chr19.block37.${pop_admx}.sim${rep}.${pcase}fun.100syn.legend \
    -L chr19.block37.${pop_admx}.sim${rep}.100fun.100syn.legend \
    -H chr19.block37.${pop_admx}.sim${rep}.all.100fun.100syn.haps.gz \
    -z

### convert the pruned haplotype file into a sparse matrix
python3 ${WD}/raresim_fix/raresim/convert.py \
    -i chr19.block37.${pop_admx}.sim${rep}.all.100fun.100syn.haps.gz \
    -o chr19.block37.${pop_admx}.sim${rep}.all.100fun.100syn.haps.sm

echo "Pruning down to ${pconf}% functional and ${pconf}% synonymous variants..."

### prune the haplotypes down to pconf % functional and pconf % synonymous variants
python3 ${WD}/raresim_fix/raresim/sim.py \
    -m chr19.block37.${pop_admx}.sim${rep}.all.100fun.100syn.haps.sm \
    --functional_bins ${WD}/mac_bin_estimates/${pop_admx}/MAC_bin_estimates_${Nsim_admx}_${pop_admx}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mac_bin_estimates/${pop_admx}/MAC_bin_estimates_${Nsim_admx}_${pop_admx}_syn_${pconf}_6bins.txt \
    --stop_threshold 10 \
    -l chr19.block37.${pop_admx}.sim${rep}.100fun.100syn.legend \
    -L chr19.block37.${pop_admx}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H chr19.block37.${pop_admx}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz \
    -z

# done