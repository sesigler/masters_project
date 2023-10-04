#!/bin/bash

# The scenarios

#################################################################
###############   Cases   # Internal Controls # Common Controls
#################################################################
#T1E          # 100% all  # 100% all	      # 100% all
#	      #		  #		      #
#	      #		  # 		      #
#################################################################
#Power	      # 120% fun  # 100% all	      #	100% all
#	      # 100% syn  #		      #
#	      #		  #		      #
#################################################################

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

cd /storage/math/projects/RAREsim/Cases/Sim_20k/${pop}/data

### replace the first lines of the Expected MAC Bins R script with the necessary variables
#Pop="Pop = '$pop'"
#Pcase="p_case = $pcase"
#Pconf="p_conf = $pconf"
#Nsim="Nsim = $nsim"
#sed -i "1c${Pop}" ${WD}/code/expected_MAC_bins_cases.R 
#sed -i "2c${Pcase}" ${WD}/code/expected_MAC_bins_cases.R 
#sed -i "3c${Pconf}" ${WD}/code/expected_MAC_bins_cases.R 
#sed -i "4c${Nsim}" ${WD}/code/expected_MAC_bins_cases.R 

### calculate the necessary expected variants (already did for NFE)
#/storage/singularity/mixtures.sif Rscript ${WD}/code/expected_MAC_bins_cases.R

### use RAREsim to prune to pfun % functional and psyn % synonymous
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
    -L ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz

# Convert -H hap file to sm
python3 ${WD}/raresim/convert.py \
    -i ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz \
    -o ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm \

# Prune fun and syn variants down again to pconf %
python3 ${WD}/raresim/sim.py \
    -m ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm \
    --functional_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_fun_${pconf}_6bins.txt \
    --synonymous_bins ${WD}/mastersProject/Input/MAC_bin_estimates_${nsim}_${pop}_syn_${pconf}_6bins.txt \
    -z \
    -l ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.${pcase}fun.${pcase}syn.legend \
    -L ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.${pconf}fun.${pconf}syn.legend \
    -H ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz

done

### replace the first lines of the Create Haplotypes R script with the necessary variables
#Ncase="Ncase = $ncase"
#sed -i "1c${Pop}" ${WD}/code/${pop}_scripts/create_haps_cases_controls_${end}.R 
#sed -i "2c${Pcase}" ${WD}/code/${pop}_scripts/create_haps_cases_controls_${end}.R 
#sed -i "3c${Pconf}" ${WD}/code/${pop}_scripts/create_haps_cases_controls_${end}.R 
#sed -i "4c${Ncase}" ${WD}/code/${pop}_scripts/create_haps_cases_controls_${end}.R 

### run the R script for case extraction and pruning back to 100% functional and common controls to  80% (all)
#/storage/singularity/mixtures.sif 
#Rscript ${WD}/code/${pop}_scripts/create_haps_cases_controls_${end}.R 