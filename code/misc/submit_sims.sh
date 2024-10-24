#!/bin/bash

# do not submit this script
# run each for loop in the terminal separately


### Subset master legend files
cd /home/math/siglersa/Sim_100k/code/
Pop=NFE

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
cd /home/math/siglersa/Sim_100k/code/
Pop=NFE
Ne=11418 #effective sample size: AFR(17469), EAS(14269), NFE(11418), SAS(14269), IAM(14269)
dl=14705483 #disease loci: AFR(14627281), EAS(14673368), NFE(14705483), SAS(14508902), IAM(14499561)

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

# submit job waiting 15 seconds in between
sbatch ./${Pop}_scripts/submit_hapgen2_${i}.sh
sleep 15

done


### Rerun duplicate seed Hapgen files
cd /home/math/siglersa/Sim_100k/code/
Pop=NFE
Ne=11418 #effective sample size: AFR(17469), EAS(14269), NFE(11418), SAS(14269), IAM(14269)
dl=14705483 #disease loci: AFR(14627281), EAS(14673368), NFE(14705483), SAS(14508902), IAM(14499561)

# Note: may have to do this more than once if duplicates still remain
# see check_sims.sh in code/misc/ on my GitHub

# make scripts for the current pop
cp ./submit_hapgen2_dups.sh ./${Pop}_scripts/submit_hapgen2_dups.sh
cp ./hapgen2_dups.sh ./${Pop}_scripts/hapgen2_dups.sh

# change "POP" to the current population
sed -i "s/POP/$Pop/g" ./${Pop}_scripts/submit_hapgen2_dups.sh 
sed -i "s/POP/$Pop/g" ./${Pop}_scripts/hapgen2_dups.sh

# change "xNEx" and "xDLx" to match the current population
sed -i "s/xNEx/$Ne/g" ./${Pop}_scripts/hapgen2_dups.sh
sed -i "s/xDLx/$dl/g" ./${Pop}_scripts/hapgen2_dups.sh

# submit job
sbatch ./${Pop}_scripts/submit_hapgen2_dups.sh


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

# Note: Check all other paramters in prune_admixed.sh are set to match the current
#       simulation scenario before running the for loop below

for i in {200..1000..100} 
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

# submit jobs
sbatch ./pruning_scripts/${Scen}/${Sub}/submit_prune_admixed_${i}.sh 

done


### Create dataset scripts
cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/code/
Scen=s2
Sub=Ncase4000_Nic4000

# Note: Check all other paramters in create_datasets_admixed.R are set to match the current
#       simulation scenario before running the for loop below

for i in {200..1000..100}
do

# make scripts for each simulation batch 
cp ./submit_create_datasets_admixed.sh  ./dataset_scripts/${Scen}/${Sub}/submit_create_datasets_admixed_${i}.sh
cp ./create_datasets_admixed.R ./dataset_scripts/${Scen}/${Sub}/create_datasets_admixed_${i}.R

# change "mysim" to the current batch number
sed -i "s/mysim/$i/g" ./dataset_scripts/${Scen}/${Sub}/submit_create_datasets_admixed_${i}.sh
sed -i "s/mysim/$i/g" ./dataset_scripts/${Scen}/${Sub}/create_datasets_admixed_${i}.R

# change "xSCENx" to the current scenario
sed -i "s/xSCENx/$Scen/g" ./dataset_scripts/${Scen}/${Sub}/submit_create_datasets_admixed_${i}.sh 
sed -i "s/xSCENx/$Scen/g" ./dataset_scripts/${Scen}/${Sub}/create_datasets_admixed_${i}.R 

# change "xSUBx" to the current sub-scenario
sed -i "s/xSUBx/$Sub/g" ./dataset_scripts/${Scen}/${Sub}/submit_create_datasets_admixed_${i}.sh 
sed -i "s/xSUBx/$Sub/g" ./dataset_scripts/${Scen}/${Sub}/create_datasets_admixed_${i}.R 

# submit job
sbatch ./dataset_scripts/${Scen}/${Sub}/submit_create_datasets_admixed_${i}.sh

done


### Create type I error scripts
cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/code/
Scen=s2
Sub=Ncase4000_Nic4000

# Note: Check all other paramters in t1e_gene_admixed.R are set to match the current
#       simulation scenario before running the for loop below

for j in {200..1000..100} 
do

# make scripts for each simulation batch
cp ./submit_t1e_gene_admixed.sh ./t1e_scripts/${Scen}/${Sub}/submit_t1e_gene_admixed_${j}.sh 
cp ./t1e_gene_admixed.R ./t1e_scripts/${Scen}/${Sub}/t1e_gene_admixed_${j}.R

# change "mysim" to the current batch number
sed -i "s/mysim/$j/g" ./t1e_scripts/${Scen}/${Sub}/submit_t1e_gene_admixed_${j}.sh
sed -i "s/mysim/$j/g" ./t1e_scripts/${Scen}/${Sub}/t1e_gene_admixed_${j}.R

# change "xSCENx" to the current scenario
sed -i "s/xSCENx/$Scen/g" ./t1e_scripts/${Scen}/${Sub}/submit_t1e_gene_admixed_${j}.sh
sed -i "s/xSCENx/$Scen/g" ./t1e_scripts/${Scen}/${Sub}/t1e_gene_admixed_${j}.R

# change "xSUBx" to the current scenario
sed -i "s/xSUBx/$Sub/g" ./t1e_scripts/${Scen}/${Sub}/submit_t1e_gene_admixed_${j}.sh
sed -i "s/xSUBx/$Sub/g" ./t1e_scripts/${Scen}/${Sub}/t1e_gene_admixed_${j}.R

# submit jobs
sbatch ./t1e_scripts/${Scen}/${Sub}/submit_t1e_gene_admixed_${j}.sh

done


### create power scripts
cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/code/
Scen=s2
Sub=Ncase4000_Nic4000

# Note: Check all other paramters in power_gene_admixed.R are set to match the current
#       simulation scenario before running the for loop below

for j in {200..1000..100} 
do

# make scripts for each simulation batch
cp ./submit_power_gene_admixed.sh ./power_scripts/${Scen}/${Sub}/submit_power_gene_admixed_${j}.sh 
cp ./power_gene_admixed.R ./power_scripts/${Scen}/${Sub}/power_gene_admixed_${j}.R

# change "mysim" to the current batch number
sed -i "s/mysim/$j/g" ./power_scripts/${Scen}/${Sub}/submit_power_gene_admixed_${j}.sh
sed -i "s/mysim/$j/g" ./power_scripts/${Scen}/${Sub}/power_gene_admixed_${j}.R

# change "xSCENx" to the current scenario
sed -i "s/xSCENx/$Scen/g" ./power_scripts/${Scen}/${Sub}/submit_power_gene_admixed_${j}.sh
sed -i "s/xSCENx/$Scen/g" ./power_scripts/${Scen}/${Sub}/power_gene_admixed_${j}.R

# change "xSUBx" to the current scenario
sed -i "s/xSUBx/$Sub/g" ./power_scripts/${Scen}/${Sub}/submit_power_gene_admixed_${j}.sh
sed -i "s/xSUBx/$Sub/g" ./power_scripts/${Scen}/${Sub}/power_gene_admixed_${j}.R

# submit jobs
sbatch ./power_scripts/${Scen}/${Sub}/submit_power_gene_admixed_${j}.sh

done


### combine type I error and power results into one file
Scen=s2
Sub=Ncase4000_Nic4000
cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/${Scen}/${Sub}/
Method=prox_ext

# keep the header from the first batch
cp ./t1e/T1e_gene_${Method}_${Scen}_${Sub}_maf0.001_100.txt ./t1e/T1e_gene_${Method}_${Scen}_${Sub}_maf0.001_NEW.txt
cp ./power/Power_gene_${Method}_${Scen}_${Sub}_maf0.001_100.txt ./power/Power_gene_${Method}_${Scen}_${Sub}_maf0.001_NEW.txt

for j in {200..1000..100} 
do

# remove the header from each file
sed -i '1d' ./t1e/T1e_gene_${Method}_${Scen}_${Sub}_maf0.001_${j}.txt
sed -i '1d' ./power/Power_gene_${Method}_${Scen}_${Sub}_maf0.001_${j}.txt

cat ./t1e/T1e_gene_${Method}_${Scen}_${Sub}_maf0.001_${j}.txt >> ./t1e/T1e_gene_${Method}_${Scen}_${Sub}_maf0.001_NEW.txt
cat ./power/Power_gene_${Method}_${Scen}_${Sub}_maf0.001_${j}.txt >> ./power/Power_gene_${Method}_${Scen}_${Sub}_maf0.001_NEW.txt

done


### calculate the type I error and power
cd /home/math/siglersa/admixed/47IAM_44NFE_5EAS_4AFR/code/
Scen=s2
Sub=Ncase4000_Nic4000
Nsim=1000

# change "xSCENx" to the current scenario
sed -i "s/xSCENx/$Scen/g" ./submit_summarize_t1e.sh
sed -i "s/xSCENx/$Scen/g" ./summarize_t1e.R

sed -i "s/xSCENx/$Scen/g" ./submit_summarize_power.sh
sed -i "s/xSCENx/$Scen/g" ./summarize_power.R

# change "xSUBx" to the current scenario
sed -i "s/xSUBx/$Sub/g" ./submit_summarize_t1e.sh
sed -i "s/xSUBx/$Sub/g" ./summarize_t1e.R

sed -i "s/xSUBx/$Sub/g" ./submit_summarize_power.sh
sed -i "s/xSUBx/$Sub/g" ./summarize_power.R

# change "xNSIMx" to the current scenario
sed -i "s/xNSIMx/$Nsim/g" ./submit_summarize_t1e.sh
sed -i "s/xNSIMx/$Nsim/g" ./submit_summarize_power.sh

# submit jobs
sbatch ./summarize_t1e.sh
sbatch ./summarize_power.sh
