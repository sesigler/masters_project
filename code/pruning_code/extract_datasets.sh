#!/bin/bash

# Extract cases, internal, and external controls from a haplotype file

pop=NFE
pruning=pruneSepRaresim
pconf=100
int_prune=100
ext_prune=80
Ncase=10000
Nint=10000
Ncc=20000

# define parameters
end=100
start=$(($end-99))

WD=/home/math/siglersa

### use RAREsim to extract the datasets from a pruned haplotype file
for rep in $(eval echo "{$start..$end}")
do

# Extract cases from pconf % pruned haplotype file
python3 ${WD}/raresim/extract.py \
    -i "${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/chr19.block37.${pop}.sim${rep}.all.${pconf}fun.${pconf}syn.haps.gz" \
    -o "${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/datasets/chr19.block37.${pop}.sim${rep}.cases.${pconf}fun.${pconf}syn.haps.gz" \
    -n ${Ncase} \
    --seed 13

# Extract internal controls from pconf % pruned haplotype REMAINDER file
python3 ${WD}/raresim/extract.py \
    -i ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/datasets/chr19.block37.${pop}.sim${rep}.cases.${pconf}fun.${pconf}syn.haps.gz-remainder \
    -o ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/datasets/chr19.block37.${pop}.sim${rep}.internal.controls.${pconf}fun.${pconf}syn.haps.gz \
    -n ${Nint} \
    --seed 13

# Extract external controls from internal controls REMAINDER file-should be all that's leftover
python3 ${WD}/raresim/extract.py \
    -i ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/datasets/chr19.block37.${pop}.sim${rep}.internal.controls.${pconf}fun.${pconf}syn.haps.gz-remainder \
    -o ${WD}/mastersProject/20K_${pop}/${pruning}/${int_prune}v${ext_prune}/datasets/chr19.block37.${pop}.sim${rep}.common.controls.${pconf}fun.${pconf}syn.haps.gz \
    -n ${Ncc} \
    --seed 13

done
