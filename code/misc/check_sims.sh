#!/bin/bash

Pop=EAS
cd /home/math/siglersa/Sim_100k/code/${Pop}_scripts/

# rm simulations.seed.txt

for i in {100..1000..100} 
do

file=hapgen2_${i}.o*

grep "seed" $file | egrep -o '[0-9.]+' >> simulations.seed.txt

done

# save duplicated seeds (507)
sort -n simulations.seed.txt | uniq -d > simulations.seed.dups.txt

# determine line number (i.e. replicate number) of duplicates (1083)
fgrep -xnf simulations.seed.dups.txt simulations.seed.txt | sort -n -t ":" -k2 > simulations.dups.lines.txt

# remove the first instance of the duplicates (don't need to be rerun) & just save the line number (576)
awk -F: 'pre==$2 { print; next }{ pre=$2 }' simulations.dups.lines.txt | cut -d ":" -f1 > simulations.rerun.lines.txt


### Rerun the replicate numbers in simulations.rerun.lines.txt & then check for duplicates again afterwards


# remove the duplicated lines from the original seed file (9425)
awk 'NR==FNR{l[$0];next;} !(FNR in l)' simulations.rerun.lines.txt simulations.seed.txt > simulations.seed.NOdups.txt

# extract the seed numbers from the rerun simulations (576)
grep "seed" hapgen2_dups.o* | egrep -o '[0-9.]+' > simulations.seed.rerun.txt

# combine the unduplicated and rerun seeds (10001)
cat simulations.seed.NOdups.txt simulations.seed.rerun.txt > simulations.seed2.txt

# double check to see if there are still duplicates (3)
sort simulations.seed2.txt | uniq -d > simulations.seed.dups2.txt

# determine line number (i.e. replicate number) of duplicates from original file (3)
fgrep -xnf simulations.seed.dups2.txt simulations.seed.txt | sort -n -t ":" -k2 > simulations.dups2.lines.txt

# just save the line number of the new duplicates
cut -d ":" -f1 simulations.dups2.lines.txt > simulations.rerun2.lines.txt


### Rerun the replicate numbers in simulations.rerun2.lines.txt & then check for duplicates again afterwards


# remove the duplicated lines from the original seed file (9998)
sort -n simulations.seed2.txt | uniq > simulations.seed2.NOdups.txt

# extract the seed numbers from the rerun simulations (576)
grep "seed" hapgen2_dups2.o* | egrep -o '[0-9.]+' > simulations.seed2.rerun.txt

# combine the unduplicated and rerun seeds (10001)
cat simulations.seed2.NOdups.txt simulations.seed2.rerun.txt > simulations.seed3.txt

# double check to see if there are still duplicates (0)
sort simulations.seed3.txt | uniq -d > simulations.seed.dups3.txt
