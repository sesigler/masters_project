############################################################################## 
# This file is used to generate the haplotype files necessary for running
# ADDITIONAL pruning scenarios for type I error calculations for proxECAT, 
# LogProx, and iECAT-O on a homogeneous population
# Note: only the additional pruned haplotype files are saved
##############################################################################

library(data.table)
library(dplyr)

source("/home/math/siglersa/mastersProject/Input/create_haps_funcs.R")
#source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/create_haps_funcs.R")


Pop = 'NFE'
p_case = 120
p_conf = 95
Nsim = 20000 

# Number of individuals in each dataset
Ncase = Nint = 5000
Ncc = 10000 #cc10k = 10000, cc5k = 5000

# Haplotype column indices
cols = 1:40000


mac_dir = '/home/math/siglersa/mastersProject/Input/'
dir_in = '/storage/math/projects/compinfo/simulations/output/20K_NFE/'
dir_out = '/home/math/siglersa/mastersProject/20K_NFE/cc10k/'

# mac_dir = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'
# dir_in = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'


### read in the expected number of functional and synonymous variants from RAREsim
exp_fun = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop, '_fun_100.txt'), header=T, sep='\t')
exp_fun_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop, '_fun_', p_conf, '.txt'), header=T, sep='\t')
exp_syn_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop, '_syn_', p_conf, '.txt'), header=T, sep='\t')


set.seed(1) # Will be different for each replicate but same for each run
#j = 1
for(j in 1:100){
  
  # read in the legend file
  leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.legend'), header=T, sep='\t')
  leg$row = 1:nrow(leg)
  
  # read in the haplotype file
  hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.100syn.haps.gz'))
  hap = as.data.frame(hap)
  
  # add allele counts to the haplotypes
  leg$count = rowSums(hap)
  
  # convert to minor allele counts
  leg$MAC = ifelse(leg$count>Nsim, 2*Nsim-leg$count, leg$count)
  
  # Select the columns for necessary data sets
  cases = sort(sample(x=cols, size = 2*Ncase, replace = FALSE))
  int = sort(sample(x=cols[! cols %in% cases], size = 2*Nint, replace = FALSE))
  cc = sort(sample(x=cols[! cols %in% c(cases, int)], size = 2*Ncc, replace = FALSE))
  
  # subset the case haplotypes (120% fun and 100% syn)
  # hap_cases = hap[, cases]
  
  # write the haplotype file for the cases (power)
  # fwrite(hap_cases, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.cases.', p_case, 'fun.100syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the legend file to the functional variants (those are the only ones we'll prune)
  leg_fun = leg %>% filter(fun=="fun")
  
  # prune the functional variants back to 100%
  rem_fun = select_var(leg_fun, exp_fun)
  hap_all_pruned = prune_var(rem_fun, hap, Nsim)
  
  # subset the pruned haplotypes for 100% pruned
  # hap_cases_pruned = hap_all_pruned[, cases]
  # hap_int = hap_all_pruned[, int]
  # hap_cc = hap_all_pruned[, cc]
  
  # write the haplotype file for the cases (type I error), internal and common controls
  # fwrite(hap_cases_pruned, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.cases.100fun.100syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  # 
  # fwrite(hap_int, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.internal.controls.100fun.100syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  # 
  # fwrite(hap_cc, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.common.controls.100fun.100syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  
  
  #### Prune back to p_conf % of the functional and synonymous variants
  
  # update the allele counts for just the common controls
  leg_cc = leg
  leg_cc$count = rowSums(hap_all_pruned)
  leg_cc$MAC = ifelse(leg_cc$count>Nsim, 2*Nsim-leg_cc$count, leg_cc$count)
  
  # subset the variants
  leg_fun_cc = leg_cc %>% filter(fun=="fun")
  leg_syn_cc = leg_cc %>% filter(fun=="syn")
  
  # prune the functional and synonymous variants of the common controls to p_conf %
  rem_cc_fun = select_var(leg_fun_cc, exp_fun_conf)
  rem_cc_syn = select_var(leg_syn_cc, exp_syn_conf)
  hap_cc_conf = prune_var(rbind(rem_cc_fun, rem_cc_syn), hap_all_pruned, Nsim)
  
  # subset the datasets for p_conf % pruning
  hap_cc_pruned = hap_cc_conf[, cc]
  # hap_case_pruned_pconf = hap_cc_conf[, cases]
  # hap_int_pruned = hap_cc_conf[, int]
  
  # write the haplotype files for p_conf % pruned
  fwrite(hap_cc_pruned, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # fwrite(hap_case_pruned_pconf, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  # 
  # fwrite(hap_int_pruned, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  
  print(j)
}
