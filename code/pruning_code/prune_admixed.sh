#!/bin/bash

# Code mostly from /storage/math/projects/compinfo/simulations/code/admixed/create_raresim_haps_2pop_v2.sh

# Simulate an admixed African American pop composed of AFR and NFE haplotypes
# Then prune haplotypes accordingly
# Note: nsim based off of max AFR and NFE individuals needed between both scenarios 1 and 2

pop1=AFR
pop2=NFE
admx_pop1=80
admx_pop2=20
Npop1=120000
Npop2=44000
Nsim_admx=82000
admx_fold=82
Nsim_hapgen=120000
haps_fold=120
pcase=160
pconf=80

# define parameters
end=100
start=$(($end-99))

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa

cd ${WD}/admixed/${admx_pop1}${pop1}_${admx_pop2}${pop2}/Sim_${admx_fold}k/${pcase}v100v${pconf}

### use RAREsim to prune to pcase % functional and 100% synonymous
for rep in $(eval echo "{$start..$end}")
do 

### Copy megan's .controls.haps.gz file into my directory then gunzip it
#cp /storage/math/projects/RAREsim/Cases/Sim_20k/${pop1}/data/chr19.block37.${pop1}.sim${rep}.controls.haps.gz .
#gunzip chr19.block37.${pop1}.sim${rep}.controls.haps.gz

# gunzip the .controls.haps.gz file
gunzip ${WD}/Sim_${haps_fold}k/${pop1}/data/chr19.block37.${pop1}.sim${rep}.${Nsim_hapgen}.controls.haps.gz

### extract the AFR haplotypes
python3 ${WD}/raresim/extract.py \
    -i ${WD}/Sim_${haps_fold}k/${pop1}/data/chr19.block37.${pop1}.sim${rep}.${Nsim_hapgen}.controls.haps \
    -o chr19.block37.${pop1}.reduced.sim${rep}.controls.haps \
    -n $Npop1 \
    --seed $rep

#gzip the gunzipped hap file
gzip ${WD}/Sim_${haps_fold}k/${pop1}/data/chr19.block37.${pop1}.sim${rep}.${Nsim_hapgen}.controls.haps

### Copy megan's .controls.haps.gz file into my directory then gunzip it
#cp /storage/math/projects/RAREsim/Cases/Sim_20k/${pop2}/data/chr19.block37.${pop2}.sim${rep}.controls.haps.gz .
#gunzip chr19.block37.${pop2}.sim${rep}.controls.haps.gz

# gunzip the .controls.haps.gz file
gunzip ${WD}/Sim_${haps_fold}k/${pop2}/data/chr19.block37.${pop2}.sim${rep}.${Nsim_hapgen}.controls.haps.gz

### extract the NFE haplotypes
python3 ${WD}/raresim/extract.py \
    -i ${WD}/Sim_${haps_fold}k/${pop2}/data/chr19.block37.${pop2}.sim${rep}.${Nsim_hapgen}.controls.haps \
    -o chr19.block37.${pop2}.reduced.sim${rep}.controls.haps \
    -n $Npop2 \
    --seed $rep

#gzip the gunzipped hap file
gzip ${WD}/Sim_${haps_fold}k/${pop2}/data/chr19.block37.${pop2}.sim${rep}.${Nsim_hapgen}.controls.haps

### combine the extracted AFR and NFE haplotype files 
paste chr19.block37.${pop1}.reduced.sim${rep}.controls.haps-sample chr19.block37.${pop2}.reduced.sim${rep}.controls.haps-sample -d " " | gzip > chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz

### convert the combined haplotype file into a sparse matrix
python3 ${WD}/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.sm

### prune the haplotypes down to pcase % functional and 100% synonymous variants
# input legend file only needs to change if admixture porportion changes
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.sm \
    --functional_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_fun_${pcase}.txt \
    --synonymous_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_syn_100.txt \
    -l ${WD}/admixed/${admx_pop1}${pop1}_${admx_pop2}${pop2}/subset_master/chr19.block37.${pop1}_${pop2}.sim${rep}.copy.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.${pcase}fun.100syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz

### extract the cases for power
#/storage/singularity/mixtures.sif python3 /storage/math/projects/RAREsim/raresim/extract.py \
#    -i chr19.block37.${pop1}-${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz \
#    -o chr19.block37.${pop1}-${pop2}.sim${rep}.cases.${pcase}fun.100syn.haps.gz \
#    -n 10000 \
#    --seed 123

### convert the pruned haplotype file into a sparse matrix
python3 ${WD}/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.sm

### prune the haplotypes down to 100% functional
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.sm \
    --f_only ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_fun_100.txt \
    -l chr19.block37.${pop1}_${pop2}.sim${rep}.${pcase}fun.100syn.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.100fun.100syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.gz

### extract the cases for type I error
#/storage/singularity/mixtures.sif python3 /storage/math/projects/RAREsim/raresim/extract.py \
#    -i chr19.block37.${pop1}-${pop2}.sim${rep}.all.100fun.100syn.haps.gz \
#    -o chr19.block37.${pop1}-${pop2}.sim${rep}.cases.100fun.100syn.haps.gz \
#    -n 10000 \
#    --seed 123

### extract the internal controls


### convert the pruned haplotype file into a sparse matrix
python3 ${WD}/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.sm

### prune the haplotypes down to pconf % functional and pconf % synonymous variants
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.sm \
    --functional_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_fun_${pconf}.txt \
    --synonymous_bins ${WD}/mac_bin_estimates/${pop1}/MAC_bin_estimates_${Nsim_admx}_${pop1}_syn_${pconf}.txt \
    -l chr19.block37.${pop1}_${pop2}.sim${rep}.100fun.100syn.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz


done

