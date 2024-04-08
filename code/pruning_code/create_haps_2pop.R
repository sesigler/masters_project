############################################################################## 
# This file is used to generate the haplotype files necessary for running
# the type I error and power calculations for proxECAT, LogProx, and iECAT-O
# on an ADMIXED population
##############################################################################
# Current set-up: Add rows of zero back in to 100% fun 100% syn and 80% fun 80% 
# syn pruned haps for power scenario
##############################################################################

library(data.table)
library(dplyr)

source("/home/math/siglersa/code/functions/create_haps_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/pruning_code/create_haps_funcs.R")

# pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
Pop1 = 'AFR'
Pop2 = 'NFE'
admx_pop1 = 80
admx_pop2 = 20
p_case = 160
p_conf = 80
Nsim = '42k'
scen = 's2'
folder = '160v100v80'
Ncase = Nic = 5000
Ncc = 10000
Nref1 = 894
Nref2 = 684
sim_params = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_', Pop1, 'ref', Nref1, '_', Pop2, 'ref', Nref2)

# Number of haplotypes in each dataset
Ncase_pop1 = Nic_pop1 = 5000 
Ncase_pop2 = Nic_pop2 = 0 

Ncc_pop1 = 16000 
Ncc_pop2 = 4000 

Nref_pop1 = Nref1*2
Nref_pop2 = Nref2*2

# Haplotype column indices
pop1_cols = 1:56000
pop2_cols = 56001:84000


dir_leg = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/pruned_haps/')
dir_in = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/pruned_haps/')
dir_out = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/', sim_params, '/', scen, '/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/')


set.seed(1) # Will be different for each replicate but same for each run
# j = 1
for(j in 1:100){
  
  # For RAREsim v2.1.1 only pruning pipeline
  leg_pcase = read.table(paste0(dir_leg, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', p_case, 'fun.100syn.legend'), header=T, sep='\t')
  leg_pcase$row = 1:nrow(leg_pcase)
  
  leg_pexp = read.table(paste0(dir_leg, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.100fun.100syn.legend'), header=T, sep='\t')
  
  leg_pconf = read.table(paste0(dir_leg, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', p_conf, 'fun.', p_conf, 'syn.legend'), header=T, sep='\t')
  
  ### For adding pruned variants back in
  hap_pcase = fread(paste0(dir_in, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.all.', p_case, 'fun.100syn.haps.gz'))
  hap_pcase = as.data.frame(hap_pcase)
  
  hap_pexp = fread(paste0(dir_in, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.all.100fun.100syn.haps.gz'))
  hap_pexp = as.data.frame(hap_pexp)
  
  hap_pconf = fread(paste0(dir_in, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.all.', p_conf, 'fun.', p_conf, 'syn.haps.gz'))
  hap_pconf = as.data.frame(hap_pconf)
  
  # Add rows of zeros back into 100% and p_conf % pruned hap files
  hap_exp_pruned = data.frame(matrix(0, nrow=nrow(leg_pcase), ncol=ncol(hap_pcase)))
  hap_exp_pruned[which(leg_pcase$id %in% leg_pexp$id),] = hap_pexp
  
  hap_all_pruned = data.frame(matrix(0, nrow=nrow(leg_pcase), ncol=ncol(hap_pcase)))
  hap_all_pruned[which(leg_pcase$id %in% leg_pconf$id),] = hap_pconf
  
  # Select the columns for necessary data sets
  cases_pop1 = sort(sample(x=pop1_cols, size = Ncase_pop1, replace = FALSE))
  ic_pop1 = sort(sample(x=pop1_cols[! pop1_cols %in% cases_pop1], size = Nic_pop1, replace = FALSE))
  cc_pop1 = sort(sample(x=pop1_cols[! pop1_cols %in% c(cases_pop1, ic_pop1)], size = Ncc_pop1, replace = FALSE))
  refs_pop1 = sort(sample(x=pop1_cols[! pop1_cols %in% c(cases_pop1, ic_pop1, cc_pop1)], size = Nref_pop1, replace = FALSE))
  
  cases_pop2 = sort(sample(x=pop2_cols, size = Ncase_pop2, replace = FALSE))
  ic_pop2 = sort(sample(x=pop2_cols[! pop2_cols %in% cases_pop2], size = Nic_pop2, replace = FALSE))
  cc_pop2 = sort(sample(x=pop2_cols[! pop2_cols %in% c(cases_pop2, ic_pop2)], size = Ncc_pop2, replace = FALSE))
  refs_pop2 = sort(sample(x=pop2_cols[! pop2_cols %in% c(cases_pop2, ic_pop2, cc_pop2)], size = Nref_pop2, replace = FALSE))
  
  # subset the case haplotypes (pcase % fun and 100% syn)
  hap_cases_pcase = hap_pcase[, c(cases_pop1, cases_pop2)]
  
  # write the haplotype file for the cases (power)
  fwrite(hap_cases_pcase, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.cases.', p_case, 'fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the pruned haplotypes for 100% pruned
  hap_cases = hap_exp_pruned[, c(cases_pop1, cases_pop2)]
  hap_int = hap_exp_pruned[, c(ic_pop1, ic_pop2)]
  hap_cc = hap_exp_pruned[, c(cc_pop1, cc_pop2)]
  hap_refs_pop1 = hap_exp_pruned[, refs_pop1]
  hap_refs_pop2 = hap_exp_pruned[, refs_pop2]
  
  # write the haplotype files for the cases, internal and common controls
  fwrite(hap_cases, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.cases.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_int, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.internal.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_cc, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.common.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_refs_pop1, paste0(dir_out, 'chr19.block37.', Pop1, '.sim', j, '.', scen, '.refs.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_refs_pop2, paste0(dir_out, 'chr19.block37.', Pop2, '.sim', j, '.', scen, '.refs.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # Subset the datasets for the p_conf % pruned haplotype
  hap_cases_pconf = hap_all_pruned[, c(cases_pop1, cases_pop2)]
  hap_int_pconf = hap_all_pruned[, c(ic_pop1, ic_pop2)]
  hap_cc_pconf = hap_all_pruned[, c(cc_pop1, cc_pop2)]
  
  # write the haplotype files for the cases, internal and common controls p_conf % pruned
  fwrite(hap_cases_pconf, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_int_pconf, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_cc_pconf, paste0(dir_out, 'chr19.block37.', Pop1, '_', Pop2, '.sim', j, '.', scen, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  print(j)
}

# Save the obs MAC bin files
# fwrite(obs_MACbin_fun_pcase, paste0(dir_out, 'obs_MACbin_fun', p_case, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_fun_exp, paste0(dir_out, 'obs_MACbin_fun', p_exp, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_fun_pconf, paste0(dir_out, 'obs_MACbin_fun', p_conf, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_syn_exp, paste0(dir_out, 'obs_MACbin_syn', p_exp, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_syn_pconf, paste0(dir_out, 'obs_MACbin_syn', p_conf, '.csv'), quote=F, row.names=F, col.names=T, sep=',')


# Vectors to store the observed variants per MAC bin
# obs_MACbin_syn_pconf <- data.frame(matrix(ncol = 9, nrow = 100))
# colnames(obs_MACbin_syn_pconf) <- c('Singletons', 'Doubletons', 'MAC.3.5', 'MAC.6.10',
#                                   'MAC.11.20', 'MAC.21.MAF0.5', 'MAF0.5.1', 'rep', 'data')
# obs_MACbin_fun_pconf = obs_MACbin_fun_exp = obs_MACbin_syn_exp = obs_MACbin_syn_pconf
# obs_MACbin_fun_pcase = obs_MACbin_fun_exp = obs_MACbin_syn_pconf = obs_MACbin_syn_exp

### read in the expected number of functional and synonymous variants from RAREsim
# exp_fun_case = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_fun_', p_case,  '.txt'), header=T, sep='\t')
# exp_syn = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_syn_', p_exp,  '.txt'), header=T, sep='\t')
# exp_fun = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_fun_', p_exp,  '.txt'), header=T, sep='\t')
# exp_fun_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_fun_', p_conf, '.txt'), header=T, sep='\t')
# exp_syn_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_syn_', p_conf, '.txt'), header=T, sep='\t')

# read in the legend file 
# leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.legend'), header=T, sep='\t')
# leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.copy.legend'), header=T, sep='\t') #For pruning OG hap file
# leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.s_only.legend'), header=T, sep='\t') # prune sequentially
# leg$row = 1:nrow(leg)

# read in the haplotype file
# hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.', p_case, 'syn.haps.gz'))
# hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.haps.gz'))
# hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.controls.haps.gz')) #For pruning OG hap file
# hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_conf, '.haps.gz')) #prune together
# hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_conf, 'fun.', p_conf, 'syn.haps.gz')) # prune separately, sequentially
# hap = as.data.frame(hap)

# add allele counts to the haplotypes
# leg$count = rowSums(hap)
# leg_pcase$count = rowSums(hap_pcase)
# leg_pcase$count_out = rowSums(hap_all_pruned)

# convert to minor allele counts
# leg$MAC = ifelse(leg$count>Nsim, 2*Nsim-leg$count, leg$count)

# subset the legend file to the functional variants (those are the only ones we'll prune)
# leg_fun = leg %>% filter(fun=="fun")
# leg_syn = leg %>% filter(fun=="syn")

# pcase % Functional
# fun120_bins1 = which(leg_fun$MAC==1)
# fun120_bins2 = which(leg_fun$MAC==2)
# fun120_bins3 = which(leg_fun$MAC>=3 & leg_fun$MAC<=5)
# fun120_bins4 = which(leg_fun$MAC>=6 & leg_fun$MAC<=exp_fun_case[4, 2])
# fun120_bins5 = which(leg_fun$MAC>=exp_fun_case[5, 1] & leg_fun$MAC<=exp_fun_case[5, 2])
# fun120_bins6 = which(leg_fun$MAC>=exp_fun_case[6, 1] & leg_fun$MAC<=exp_fun_case[6, 2])
# fun120_bins7 = which(leg_fun$MAC>=exp_fun_case[7, 1] & leg_fun$MAC<=exp_fun_case[7, 2])
# 
# funBins_pcase = c(length(fun120_bins1), length(fun120_bins2),length(fun120_bins3),
#                   length(fun120_bins4), length(fun120_bins5), length(fun120_bins6),
#                   length(fun120_bins7))
# 
# obs_MACbin_fun_pcase[j, 1:7] <- funBins_pcase
# obs_MACbin_fun_pcase[j, 8] <- j
# obs_MACbin_fun_pcase[j, 9] <- paste0('RAREsim functional-', p_case, '%')

# MAC_ests_fun_120 = exp_fun_case
# MAC_ests_fun_120$Observed = c(length(fun120_bins1), length(fun120_bins2),length(fun120_bins3),
#                               length(fun120_bins4), length(fun120_bins5), length(fun120_bins6),
#                               length(fun120_bins7))
# fwrite(MAC_ests_fun_120, paste0(dir_out, 'MAC_bin_ests_sim', j, '_fun_', p_case, '_', Ncc, '.csv'),
#        quote=F, row.names=F, col.names=T, sep=',')

# 100 % Synonymous
# syn100_bins1 = which(leg_syn$MAC==1)
# syn100_bins2 = which(leg_syn$MAC==2)
# syn100_bins3 = which(leg_syn$MAC>=3 & leg_syn$MAC<=5)
# syn100_bins4 = which(leg_syn$MAC>=6 & leg_syn$MAC<=exp_syn[4, 2])
# syn100_bins5 = which(leg_syn$MAC>=exp_syn[5, 1] & leg_syn$MAC<=exp_syn[5, 2])
# syn100_bins6 = which(leg_syn$MAC>=exp_syn[6, 1] & leg_syn$MAC<=exp_syn[6, 2])
# syn100_bins7 = which(leg_syn$MAC>=exp_syn[7, 1] & leg_syn$MAC<=exp_syn[7, 2])
# 
# synBins_exp = c(length(syn100_bins1), length(syn100_bins2),length(syn100_bins3),
#                 length(syn100_bins4), length(syn100_bins5), length(syn100_bins6),
#                 length(syn100_bins7))
# 
# obs_MACbin_syn_exp[j, 1:7] <- synBins_exp
# obs_MACbin_syn_exp[j, 8] <- j
# obs_MACbin_syn_exp[j, 9] <- paste0('RAREsim synonymous-', p_exp, '%')

# 100 % Functional
# fun_bins1 = which(leg_fun$MAC==1)
# fun_bins2 = which(leg_fun$MAC==2)
# fun_bins3 = which(leg_fun$MAC>=3 & leg_fun$MAC<=5)
# fun_bins4 = which(leg_fun$MAC>=6 & leg_fun$MAC<=exp_fun[4, 2])
# fun_bins5 = which(leg_fun$MAC>=exp_fun[5, 1] & leg_fun$MAC<=exp_fun[5, 2])
# fun_bins6 = which(leg_fun$MAC>=exp_fun[6, 1] & leg_fun$MAC<=exp_fun[6, 2])
# fun_bins7 = which(leg_fun$MAC>=exp_fun[7, 1] & leg_fun$MAC<=exp_fun[7, 2])
# 
# funBins_pexp = c(length(fun_bins1), length(fun_bins2),length(fun_bins3),
#                   length(fun_bins4), length(fun_bins5), length(fun_bins6),
#                   length(fun_bins7))
# 
# obs_MACbin_fun_exp[j, 1:7] <- funBins_pexp
# obs_MACbin_fun_exp[j, 8] <- j
# obs_MACbin_fun_exp[j, 9] <- paste0('RAREsim functional-', p_exp, '%')

# MAC_ests_syn_100 = exp_syn
# MAC_ests_syn_100$Observed = c(length(syn100_bins1), length(syn100_bins2),length(syn100_bins3),
#                               length(syn100_bins4), length(syn100_bins5), length(syn100_bins6),
#                               length(syn100_bins7))
# fwrite(MAC_ests_syn_100, paste0(dir_out, 'MAC_bin_ests_sim', j, '_syn_100_', Ncc, '.csv'),
#        quote=F, row.names=F, col.names=T, sep=',')


# prune the functional variants back to 100%
# rem_fun = select_var(leg_fun, exp_fun)
# hap_all_pruned = prune_var(rem_fun, hap, Nsim)

# prune the synonymous variants to 100%
# rem_syn = select_var(leg_syn, exp_syn)
# hap_all_pruned = prune_var(rem_syn, hap, Nsim)

# PRUNE OG Hap file to 100% fun and 100% syn
# rem_fun = select_var(leg_fun, exp_fun)
# rem_syn = select_var(leg_syn, exp_syn)
# hap_exp = prune_var(rbind(rem_fun, rem_syn), hap, Nsim)

# hap_cases = hap_pcase[, cases]
# hap_int = hap_pcase[, int]
# hap_cc = hap_pcase[, cc]
# hap_cases = hap_exp[, cases]
# hap_int = hap_exp[, int]
# hap_cc = hap_exp[, cc]

# hap_all_pruned = hap # just need if going from 100% to 80%
# hap_all_pruned = hap_exp

# pconf % Functional
# fun_bins1 = which(leg_fun$MAC==1)
# fun_bins2 = which(leg_fun$MAC==2)
# fun_bins3 = which(leg_fun$MAC>=3 & leg_fun$MAC<=5)
# fun_bins4 = which(leg_fun$MAC>=6 & leg_fun$MAC<=exp_fun_conf[4, 2])
# fun_bins5 = which(leg_fun$MAC>=exp_fun_conf[5, 1] & leg_fun$MAC<=exp_fun_conf[5, 2])
# fun_bins6 = which(leg_fun$MAC>=exp_fun_conf[6, 1] & leg_fun$MAC<=exp_fun_conf[6, 2])
# fun_bins7 = which(leg_fun$MAC>=exp_fun_conf[7, 1] & leg_fun$MAC<=exp_fun_conf[7, 2])
# 
# funBins_pconf = c(length(fun_bins1), length(fun_bins2),length(fun_bins3),
#                   length(fun_bins4), length(fun_bins5), length(fun_bins6),
#                   length(fun_bins7))
# 
# obs_MACbin_fun_pconf[j, 1:7] <- funBins_pconf
# obs_MACbin_fun_pconf[j, 8] <- j
# obs_MACbin_fun_pconf[j, 9] <- paste0('RAREsim functional-', p_conf, '%')

# pconf % Synonymous
# syn_bins1 = which(leg_syn$MAC==1)
# syn_bins2 = which(leg_syn$MAC==2)
# syn_bins3 = which(leg_syn$MAC>=3 & leg_syn$MAC<=5)
# syn_bins4 = which(leg_syn$MAC>=6 & leg_syn$MAC<=exp_syn_conf[4, 2])
# syn_bins5 = which(leg_syn$MAC>=exp_syn_conf[5, 1] & leg_syn$MAC<=exp_syn_conf[5, 2])
# syn_bins6 = which(leg_syn$MAC>=exp_syn_conf[6, 1] & leg_syn$MAC<=exp_syn_conf[6, 2])
# syn_bins7 = which(leg_syn$MAC>=exp_syn_conf[7, 1] & leg_syn$MAC<=exp_syn_conf[7, 2])
# 
# synBins_pconf = c(length(syn_bins1), length(syn_bins2),length(syn_bins3),
#                   length(syn_bins4), length(syn_bins5), length(syn_bins6),
#                   length(syn_bins7))
# 
# obs_MACbin_syn_pconf[j, 1:7] <- synBins_pconf
# obs_MACbin_syn_pconf[j, 8] <- j
# obs_MACbin_syn_pconf[j, 9] <- paste0('RAREsim synonymous-', p_conf, '%')

# 100% Functional
# ap_fun = hap_all_pruned[leg_fun$row, ]
# ap_fun = hap_exp[leg_fun$row, ]
# ap_fun_sums = rowSums(ap_fun)
# 
# fun100_bins1 = which(ap_fun_sums==1)
# fun100_bins2 = which(ap_fun_sums==2)
# fun100_bins3 = which(ap_fun_sums>=3 & ap_fun_sums<=5)
# fun100_bins4 = which(ap_fun_sums>=6 & ap_fun_sums<=exp_fun[4, 2])
# fun100_bins5 = which(ap_fun_sums>=exp_fun[5, 1] & ap_fun_sums<=exp_fun[5, 2])
# fun100_bins6 = which(ap_fun_sums>=exp_fun[6, 1] & ap_fun_sums<=exp_fun[6, 2])
# fun100_bins7 = which(ap_fun_sums>=exp_fun[7, 1] & ap_fun_sums<=exp_fun[7, 2])
# 
# funBins_exp = c(length(fun100_bins1), length(fun100_bins2),length(fun100_bins3),
#                 length(fun100_bins4), length(fun100_bins5), length(fun100_bins6),
#                 length(fun100_bins7))
# 
# obs_MACbin_fun_exp[j, 1:7] <- funBins_exp
# obs_MACbin_fun_exp[j, 8] <- j
# obs_MACbin_fun_exp[j, 9] <- paste0('RAREsim functional-100%')

# MAC_ests_fun_100 = exp_fun
# MAC_ests_fun_100$Observed = c(length(fun100_bins1), length(fun100_bins2),length(fun100_bins3),
#                               length(fun100_bins4), length(fun100_bins5), length(fun100_bins6),
#                               length(fun100_bins7))
# fwrite(MAC_ests_fun_100, paste0(dir_out, 'MAC_bin_ests_sim', j, '_fun_100_', Ncc, '.csv'),
#        quote=F, row.names=F, col.names=T, sep=',')

# 100% Synonymous
# ap_syn = hap_exp[leg_syn$row, ]
# ap_syn_sums = rowSums(ap_syn)
# 
# syn100_bins1 = which(ap_syn_sums==1)
# syn100_bins2 = which(ap_syn_sums==2)
# syn100_bins3 = which(ap_syn_sums>=3 & ap_syn_sums<=5)
# syn100_bins4 = which(ap_syn_sums>=6 & ap_syn_sums<=exp_syn[4, 2])
# syn100_bins5 = which(ap_syn_sums>=exp_syn[5, 1] & ap_syn_sums<=exp_syn[5, 2])
# syn100_bins6 = which(ap_syn_sums>=exp_syn[6, 1] & ap_syn_sums<=exp_syn[6, 2])
# syn100_bins7 = which(ap_syn_sums>=exp_syn[7, 1] & ap_syn_sums<=exp_syn[7, 2])
# 
# synBins_exp = c(length(syn100_bins1), length(syn100_bins2),length(syn100_bins3),
#                 length(syn100_bins4), length(syn100_bins5), length(syn100_bins6),
#                 length(syn100_bins7))
# 
# obs_MACbin_syn_exp[j, 1:7] <- synBins_exp
# obs_MACbin_syn_exp[j, 8] <- j
# obs_MACbin_syn_exp[j, 9] <- paste0('RAREsim synonymous-100%')

#### Prune back to p_conf% of the functional and synonymous variants

# update the allele counts for just the common controls
# leg_cc = leg
# leg_cc$count = rowSums(hap_all_pruned)
# leg_cc$MAC = ifelse(leg_cc$count>Nsim, 2*Nsim-leg_cc$count, leg_cc$count)

# subset the variants
# leg_fun_cc = leg_cc %>% filter(fun=="fun")
# leg_syn_cc = leg_cc %>% filter(fun=="syn")

# prune the functional and synonymous variants of the common controls to p_conf%
# rem_cc_fun = select_var(leg_fun_cc, exp_fun_conf)
# rem_cc_syn = select_var(leg_syn_cc, exp_syn_conf)
# hap_cc_conf = prune_var(rbind(rem_cc_fun, rem_cc_syn), hap_all_pruned, Nsim)

# subset the datasets for p_conf% pruning
# hap_cc_pruned = hap_cc_conf[, cc]
# hap_case_pconf = hap_cc_conf[, cases]
# hap_int_pruned = hap_cc_conf[, int]

# write the haplotype files for 80% pruned
# fwrite(hap_cc_pruned, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
#        quote=F, row.names=F, col.names=F, sep=' ')
# 
# fwrite(hap_case_pconf, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
#        quote=F, row.names=F, col.names=F, sep=' ')
# 
# fwrite(hap_int_pruned, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
#        quote=F, row.names=F, col.names=F, sep=' ')

# p_conf % Synonymous
# hap_cc_syn = hap_cc_conf[leg_syn_cc$row, ]
# hap_syn_sums = rowSums(hap_cc_syn)
# 
# syn80_bins1 = which(hap_syn_sums==1)
# syn80_bins2 = which(hap_syn_sums==2)
# syn80_bins3 = which(hap_syn_sums>=3 & hap_syn_sums<=5)
# syn80_bins4 = which(hap_syn_sums>=6 & hap_syn_sums<=exp_syn_conf[4, 2])
# syn80_bins5 = which(hap_syn_sums>=exp_syn_conf[5, 1] & hap_syn_sums<=exp_syn_conf[5, 2])
# syn80_bins6 = which(hap_syn_sums>=exp_syn_conf[6, 1] & hap_syn_sums<=exp_syn_conf[6, 2])
# syn80_bins7 = which(hap_syn_sums>=exp_syn_conf[7, 1] & hap_syn_sums<=exp_syn_conf[7, 2])
# 
# synBins_pconf = c(length(syn80_bins1), length(syn80_bins2),length(syn80_bins3),
#                   length(syn80_bins4), length(syn80_bins5), length(syn80_bins6),
#                   length(syn80_bins7))
# 
# obs_MACbin_syn_pconf[j, 1:7] <- synBins_pconf
# obs_MACbin_syn_pconf[j, 8] <- j
# obs_MACbin_syn_pconf[j, 9] <- paste0('RAREsim synonymous-', p_conf, '%')

# MAC_ests_syn_80 = exp_syn_conf
# MAC_ests_syn_80$Observed = c(length(syn80_bins1), length(syn80_bins2),length(syn80_bins3),
#                              length(syn80_bins4), length(syn80_bins5), length(syn80_bins6),
#                              length(syn80_bins7))
# fwrite(MAC_ests_syn_80, paste0(dir_out, 'MAC_bin_ests_sim', j, '_syn_', p_conf, '_', Ncc, '.csv'),
#        quote=F, row.names=F, col.names=T, sep=',')

# p_conf % Functional
# hap_cc_fun = hap_cc_conf[leg_fun_cc$row, ]
# hap_fun_sums = rowSums(hap_cc_fun)
# 
# fun80_bins1 = which(hap_fun_sums==1)
# fun80_bins2 = which(hap_fun_sums==2)
# fun80_bins3 = which(hap_fun_sums>=3 & hap_fun_sums<=5)
# fun80_bins4 = which(hap_fun_sums>=6 & hap_fun_sums<=exp_fun_conf[4, 2])
# fun80_bins5 = which(hap_fun_sums>=exp_fun_conf[5, 1] & hap_fun_sums<=exp_fun_conf[5, 2])
# fun80_bins6 = which(hap_fun_sums>=exp_fun_conf[6, 1] & hap_fun_sums<=exp_fun_conf[6, 2])
# fun80_bins7 = which(hap_fun_sums>=exp_fun_conf[7, 1] & hap_fun_sums<=exp_fun_conf[7, 2])
# 
# funBins_pconf = c(length(fun80_bins1), length(fun80_bins2),length(fun80_bins3),
#                   length(fun80_bins4), length(fun80_bins5), length(fun80_bins6),
#                   length(fun80_bins7))
# 
# obs_MACbin_fun_pconf[j, 1:7] <- funBins_pconf
# obs_MACbin_fun_pconf[j, 8] <- j
# obs_MACbin_fun_pconf[j, 9] <- paste0('RAREsim functional-', p_conf, '%')

# MAC_ests_fun_80 = exp_fun_conf
# MAC_ests_fun_80$Observed = c(length(fun80_bins1), length(fun80_bins2),length(fun80_bins3),
#                              length(fun80_bins4), length(fun80_bins5), length(fun80_bins6),
#                              length(fun80_bins7))
# fwrite(MAC_ests_fun_80, paste0(dir_out, 'MAC_bin_ests_sim', j, '_fun_', p_conf, '_', Ncc, '.csv'),
#        quote=F, row.names=F, col.names=T, sep=',')
