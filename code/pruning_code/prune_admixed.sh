#!/bin/bash

# Code mostly from /storage/math/projects/compinfo/simulations/code/admixed/create_raresim_haps_2pop_v2.sh

# Simulate an admixed African American pop composed of AFR and NFE haplotypes
# Then prune haplotypes accordingly

pop1=AFR
pop2=NFE
nsim=23000
pcase=160
pconf=80

# define parameters
#end=100
#start=$(($end-99))

#WD=/storage/math/projects/compinfo/simulations
WD=/home/math/siglersa

### replace the first lines of the Expected MAC Bins R script with the necessary variables
#Pop="Pop = '$pop1'"
#Pcase="p_case = $pcase"
#Pconf="p_conf = $pconf"
#Nsim="Nsim = $nsim"
#sed -i "1c${Pop}" ${WD}/code/expected_MAC_bins_cases.R 
#sed -i "2c${Pcase}" ${WD}/code/expected_MAC_bins_cases.R 
#sed -i "3c${Pconf}" ${WD}/code/expected_MAC_bins_cases.R 
#sed -i "4c${Nsim}" ${WD}/code/expected_MAC_bins_cases.R 

### calculate the necessary expected variants (just need AFR)
#/storage/singularity/mixtures.sif Rscript ${WD}/code/expected_MAC_bins_cases.R


### use RAREsim to prune to pcase % functional and 100% synonymous
#for rep in $(eval echo "{$start..$end}")
for rep in {1..10}
do 

### extract the AFR haplotypes
cd /storage/math/projects/RAREsim/Cases/Sim_20k/${pop1}/data
python3 ${WD}/raresim/extract.py \
    -i chr19.block37.${pop1}.sim${rep}.controls.haps.gz \
    -o ${WD}/admixed/${pop1}_${pop2}_pops/chr19.block37.${pop1}.reduced.sim${rep}.controls.haps \
    -n 37000 \
    --seed 123

### extract the NFE haplotypes
cd /storage/math/projects/RAREsim/Cases/Sim_20k/${pop2}/data
python3 ${WD}/raresim/extract.py \
    -i chr19.block37.${pop2}.sim${rep}.controls.haps.gz \
    -o ${WD}/admixed/${pop1}_${pop2}_pops/chr19.block37.${pop2}.reduced.sim${rep}.controls.haps \
    -n 9000 \
    --seed 123

### combine the extracted AFR and NFE haplotype files 
cd ${WD}/admixed/${pop1}_${pop2}_pops
paste chr19.block37.${pop1}.reduced.sim${rep}.controls.haps chr19.block37.${pop2}.reduced.sim${rep}.controls.haps -d " " | gzip > chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz

### convert the combined haplotype file into a sparse matrix
python3 ${WD}/raresim/convert.py \
    -i chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz \
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz.sm

### prune the haplotypes down to pcase % functional and 100% synonymous variants
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.controls.haps.gz.sm \
    --functional_bins ${WD}/mastersProject/input/MAC_bin_estimates_${nsim}_${pop1}_fun_${pcase}.txt \
    --synonymous_bins ${WD}/mastersProject/input/MAC_bin_estimates_${nsim}_${pop1}_syn_100.txt \
    -l chr19.block37.${pop1}_${pop2}.sim${rep}.copy.legend \
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
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz.sm

### prune the haplotypes down to 100% functional
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pcase}fun.100syn.haps.gz.sm \
    --f_only ${WD}/mastersProject/input/MAC_bin_estimates_${nsim}_${pop1}_fun_100.txt \
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
    -o chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.gz.sm

### prune the haplotypes down to pconf % functional and pconf % synonymous variants
python3 ${WD}/raresim/sim.py \
    -m chr19.block37.${pop1}_${pop2}.sim${rep}.all.100fun.100syn.haps.gz.sm \
    --functional_bins ${WD}/mastersProject/input/MAC_bin_estimates_${nsim}_${pop1}_fun_${pconf}.txt \
    --synonymous_bins ${WD}/mastersProject/input/MAC_bin_estimates_${nsim}_${pop1}_syn_${pconf}.txt \
    -l chr19.block37.${pop1}_${pop2}.sim${rep}.100fun.100syn.legend \
    -L chr19.block37.${pop1}_${pop2}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H chr19.block37.${pop1}_${pop2}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz


done

