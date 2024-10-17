#!/bin/bash

# do not submit this script
# run each for loop in the terminal separately

cd /home/math/siglersa/Sim_100k/code/
Pop=IAM
Ne=14269 #effective sample size: AFR(17469), EAS(14269), NFE(11418), SAS(14269), IAM(14269)
dl=14499561 #disease loci: AFR(14627281), EAS(14673368), NFE(14705483), SAS(14508902), IAM(14499561)

### Subset master legend files
for i in {200..1000..100}
do

# make scripts for each simulation batch 
cp ./submit_subset_master.sh  ./${Pop}_scripts/submit_subset_master_${i}.sh
cp ./subset_master.R ./${Pop}_scripts/subset_master_${i}.R

# change "mysim" to the current batch number
sed -i "s/mysim/$i/g" ./${Pop}_scripts/submit_subset_master_${i}.sh
sed -i "s/mysim/$i/g" ./${Pop}_scripts/subset_master_${i}.R

# change "POP" to the current population
sed -i "s/POP/$Pop/g" ./${Pop}_scripts/submit_subset_master_${i}.sh 
sed -i "s/POP/$Pop/g" ./${Pop}_scripts/subset_master_${i}.R

# submit job
sbatch ./${Pop}_scripts/submit_subset_master_${i}.sh

done


### Create HAPGEN scripts
for i in {200..1000..100}
do

# make scripts for each simulation batch
cp ./submit_hapgen2.sh ./${Pop}_scripts/submit_hapgen2_${i}.sh
cp ./hapgen2.sh ./${Pop}_scripts/hapgen2_${i}.sh

# change "mysim" to the current batch number
sed -i "s/mysim/$i/g" ./${Pop}_scripts/submit_hapgen2_${i}.sh
sed -i "s/mysim/$i/g" ./${Pop}_scripts/hapgen2_${i}.sh

# change "POP" to the current population
sed -i "s/POP/$Pop/g" ./${Pop}_scripts/submit_hapgen2_${i}.sh 
sed -i "s/POP/$Pop/g" ./${Pop}_scripts/hapgen2_${i}.sh

# change "xNEx" and "xDLx" to match the current population
sed -i "s/xNEx/$Ne/g" ./${Pop}_scripts/hapgen2_${i}.sh
sed -i "s/xDLx/$dl/g" ./${Pop}_scripts/hapgen2_${i}.sh

# submit jobw aiting 15 seconds in between
sbatch ./${Pop}_scripts/submit_hapgen2_${i}.sh
sleep 15

done

### Subset admixed master legend files

cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/code/

for i in {200..1000..100}
do

# make scripts for each simulation batch 
cp ./submit_subset_master_admixed.sh  ./master_admixed_scripts/submit_subset_master_admixed_${i}.sh
cp ./subset_master_admixed.R ./master_admixed_scripts/subset_master_admixed_${i}.R

# change "mysim" to the current batch number
sed -i "s/mysim/$i/g" ./master_admixed_scripts/submit_subset_master_admixed_${i}.sh
sed -i "s/mysim/$i/g" ./master_admixed_scripts/subset_master_admixed_${i}.R

# submit job
sbatch ./master_admixed_scripts/submit_subset_master_admixed_${i}.sh

done

### Create RAREsim scripts
cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/code/
Scen=s2
Sub=Ncase4000_Nic4000

for i in {100..1000..100} 
do

# make scripts for each simulation batch
cp ./submit_prune_admixed.sh ./pruning_scripts/${Scen}/${Sub}/submit_prune_admixed_${i}.sh 
cp ./prune_admixed.sh ./pruning_scripts/${Scen}/${Sub}/prune_admixed_${i}.sh

# change "mysim" to the current batch number
sed -i "s/mysim/$i/g" ./pruning_scripts/${Scen}/${Sub}/submit_prune_admixed_${i}.sh 
sed -i "s/mysim/$i/g" ./pruning_scripts/${Scen}/${Sub}/prune_admixed_${i}.sh 

# change "xSCENx" to the current scenario
sed -i "s/xSCENx/$Scen/g" ./pruning_scripts/${Scen}/${Sub}/submit_prune_admixed_${i}.sh 
sed -i "s/xSCENx/$Scen/g" ./pruning_scripts/${Scen}/${Sub}/prune_admixed_${i}.sh 

# change "xSUBx" to the current sub-scenario
sed -i "s/xSUBx/$Sub/g" ./pruning_scripts/${Scen}/${Sub}/submit_prune_admixed_${i}.sh 
sed -i "s/xSUBx/$Sub/g" ./pruning_scripts/${Scen}/${Sub}/prune_admixed_${i}.sh 

# submit jobs waiting 15 seconds in between
sbatch ./pruning_scripts/${Scen}/${Sub}/submit_prune_admixed_${i}.sh 
#sleep 15

done


### create type I error and power scripts
for j in {100..1000..100} 
do

# make scripts for each simulation batch
cp ./methods.R ./${Pop}_scripts/methods_${j}.R
cp ./methods.sh ./${Pop}_scripts/methods_${j}.sh 

# change "mysim" to the current batch number
sed -i "s/mysim/$j/g" ./${Pop}_scripts/methods_${j}.R
sed -i "s/mysim/$j/g" ./${Pop}_scripts/methods_${j}.sh

# change "pop" to the current population
sed -i "s/pop/$Pop/g" ./${Pop}_scripts/methods_${j}.R
sed -i "s/pop/$Pop/g" ./${Pop}_scripts/methods_${j}.sh

sbatch ./${Pop}_scripts/methods_${j}.sh

done


### combine type I error and power results into one file
cd /storage/math/projects/compinfo/simulations/results/

# keep the header from the first batch
cp ./20K_${Pop}/Type1_${Pop}_maf0.001_100.txt ./20K_${Pop}/Type1_${Pop}_maf0.001_NEW.txt
cp ./20K_${Pop}/Power_${Pop}_maf0.001_100.txt ./20K_${Pop}/Power_${Pop}_maf0.001_NEW.txt

for j in {200..1000..100} 
do

# remove the header from each file
sed -i '1d' ./20K_${Pop}/Type1_${Pop}_maf0.001_${j}.txt
sed -i '1d' ./20K_${Pop}/Power_${Pop}_maf0.001_${j}.txt

cat ./20K_${Pop}/Type1_${Pop}_maf0.001_${j}.txt >> ./20K_${Pop}/Type1_${Pop}_maf0.001_NEW.txt
cat ./20K_${Pop}/Power_${Pop}_maf0.001_${j}.txt >> ./20K_${Pop}/Power_${Pop}_maf0.001_NEW.txt

done


### calculate the type I error and power
cd /storage/math/projects/compinfo/simulations/code/

# make scripts for each population
cp ./summarize_results.sh ./${Pop}_scripts/summarize_results.sh 
cp ./summarize_results.R ./${Pop}_scripts/summarize_results.R

# change "pop" to the current population
sed -i "s/pop/$Pop/g" ./${Pop}_scripts/summarize_results.sh 
sed -i "s/pop/$Pop/g" ./${Pop}_scripts/summarize_results.R

sbatch ./${Pop}_scripts/summarize_results.sh
