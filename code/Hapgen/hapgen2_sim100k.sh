#!/bin/bash

end=100
start=$(($end-99))

### variables to be updated per simulation
pop=IAM #population: AFR(African), EAS(Eastern Asian), NFE(Non-Finnish European), SAS(South Asian), IAM(Indigenous American)
NE=14269 #effective sample size: AFR(17469), EAS(14269), NFE(11418), SAS(14269)
Nsim=100000 #number of individuals
num=19 #chromosome number
b=37 #block number

### input variables
dir=/home/math/siglersa
WD=/storage/math/projects/compinfo/simulations
Hap=${dir}/IAM_haps/master_legend/IAM_chr19_block37_coding_region.hap
Map=${WD}/input/genetic_map_chr${num}_combined.txt
#DL=$(cat ${WD}/input/${pop}_blocks/${pop}block_disease_loci.txt | sed -n "38p") #38p for NFE
### Calculated Disease locus by hand (random non-zero position)

for i in $(eval echo "{$start..$end}")
do

### output variables
# Leg is from subset_master.R
Leg=${dir}/Sim_100k/IAM/subset_master/chr${num}.block${b}.${pop}.sim${i}.legend
Output=${dir}/Sim_100k/${pop}/data/chr${num}.block${b}.${pop}.sim${i}

### save a copy of the legend file
cp $Leg $Output.copy.legend

### simulate with HAPGEN2
#/storage/singularity/mixtures.sif 
hapgen2 -h $Hap \
-m $Map \
-l $Leg \
-o $Output.gz \
-n $Nsim 0 \
-dl 14884751 1 1.0 1.0 \
-Ne $NE \
-no_gens_output

### remove unnecessary files produced by Hapgen
rm $Output.cases.*

done