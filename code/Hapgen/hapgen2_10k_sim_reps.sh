#!/bin/bash

end=mysim
start=$(($end-99))

### variables to be updated per simulation
pop=AFR #population: AFR(African), EAS(Eastern Asian), NFE(Non-Finnish European), SAS(South Asian)
NE=17469 #effective sample size: AFR(17469), EAS(14269), NFE(11418), SAS(14269)
Nsim=120000 #number of individuals
Nind=120 #short version of Nsim for folder
num=19 #chromosome number
b=37 #block number

### input variables
MD=/home/math/siglersa
WD=/storage/math/projects/compinfo/simulations
Hap=${WD}/input/${pop}_blocks/${pop}_Block${b}_CDS_ref_added.hap
Map=${WD}/input/genetic_map_chr${num}_combined.txt
#Master=${WD}/input/chr${num}.block${b}.${pop}.master.legend
#DL=$(cat ${WD}/input/${pop}_blocks/${pop}block_disease_loci.txt | sed -n "38p") #38p for NFE
### Calculated Disease locus by hand (random non-zero position)

for i in $(eval echo "{$start..$end}")
do

### output variables
Leg=/storage/math/projects/RAREsim/Cases/Sim_20k/${pop}/data/chr${num}.block${b}.${pop}.sim${i}.copy.legend
Output=${MD}/Sim_${Nind}k/${pop}/data/chr${num}.block${b}.${pop}.sim${i}.${Nsim}
#Output=${WD}/output/multiple_cohorts/chr${num}.block${b}.${pop}.sim${i}.${Nsim}

### save a copy of the legend file
cp $Leg $Output.copy.legend

### simulate with HAPGEN2
#/storage/singularity/mixtures.sif 
hapgen2 -h $Hap \
-m $Map \
-l $Leg \
-o $Output.gz \
-n $Nsim 0 \
-dl 14627281 1 1.0 1.0 \
-Ne $NE \
-no_gens_output

### remove unnecessary files produced by Hapgen
rm $Output.cases.*

done