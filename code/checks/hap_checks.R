### Check why extract.py not working

library(data.table)
library(dplyr)

# source("/home/math/siglersa/mastersProject/Input/create_haps_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/create_haps_funcs.R")

Pop = 'NFE'
p_case = 100
# p_exp = 100
p_conf = 80
Nsim = 20000 
pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim
folder = '100v80'
int_prune = 100
ext_prune = 100

# Number of individuals in each dataset
Ncase = Nint = 5000
Ncc = 10000 #cc10k = 10000, cc5k = 5000

# Haplotype column indices
cols = 1:40000

# mac_dir = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/input/'
dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/output/'
j=1
hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.', p_case, 'syn.haps.gz'))
hap = as.data.frame(hap)

hap_pconf = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_conf, 'fun.', p_conf, 'syn.haps.gz'))
hap_pconf = as.data.frame(hap_pconf)

hap_examp = fread(paste0(dir_in, 'Simulated_80k_9.controls.haps'))
hap_examp_gz = fread(paste0(dir_in, 'Simulated_80k_9.controls.haps.gz'))
