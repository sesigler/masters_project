#!/bin/bash

STARTTIME=$(date +%s)


### variables to be updated per simulation
pop=AFR #population: AFR(African), EAS(Eastern Asian), NFE(Non-Finnish European), SAS(South Asian)
NE=17469 #effective sample size: AFR(17469), EAS(14269), NFE(11418), SAS(14269)
Nsim=20000 #number of individuals
num=19 #chromosome number
b=37 #block number

### input variables
WD_input=/storage/math/projects/compinfo/simulations
WD=/storage/math/projects/RAREsim/Cases/Sim_20k/AFR
Hap=${WD_input}/input/${pop}_blocks/${pop}_Block${b}_CDS_ref_added.hap
Map=${WD_input}/input/genetic_map_chr${num}_combined.txt
Master=${WD_input}/input/chr${num}.block${b}.${pop}.master.legend

# DL=$(cat ${WD_input}/input/${pop}_blocks/${pop}block_disease_loci.txt | sed -n "38p") 
### Calculated Disease locus by hand (random non-zero position)

### For rarest:
pyWD=/storage/math/projects/compinfo/simulations/code

#TIME7=$(date +%s)

### subset single legend files from the master legend file

Rscript ${WD}/code/subset_master.R

#TIME8=$(date +%s)
#echo "It took $(($TIME8 - $TIME7)) seconds to subset 100 single legend files from the master legend file"

for i in {2..1000}
do

STARTTIME2=$(date +%s)

### output variables
Leg=${WD}/data/chr${num}.block${b}.${pop}.sim${i}.legend
Output=${WD}/data/chr${num}.block${b}.${pop}.sim${i}

## The extra columns get erased  when  we run through HAPGEN2  -  make a copy
#cp ${WD}/data/chr${num}.block${b}.${pop}.sim${i}.legend ${WD}/data/chr${num}.block${b}.${pop}.sim${i}.copy.legend

### simulate with HAPGEN2
#hapgen2 -h $Hap \
#-m $Map \
#-l $Leg \
#-o $Output.gz \
#-n $Nsim 0 \
#-dl 14627281 1 1.0 1.0 \
#-no_gens_output \
#-Ne $NE

#TIME1=$(date +%s)
#echo "It took $(($TIME1 - $STARTTIME2)) seconds to run Hapgen2 for simulation $i"


### And create the sparse matrix:
#/storage/singularity/mixtures.sif 
python3 ${pyWD}/convert.py \
    -i ${WD}/data/chr19.block37.${pop}.sim${i}.controls.haps.gz \
    -o ${WD}/data/chr19.block37.${pop}.sim${i}.controls.haps.sm


done

rm ${WD}/data/*.cases.*
rm ${WD}/data/*.sample
