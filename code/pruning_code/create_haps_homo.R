############################################################################## 
# This file is used to generate the haplotype files necessary for running
# the type I error and power calculations for proxECAT, LogProx, and iECAT-O
# on a homogeneous population
##############################################################################
# Current set-up: Subset the datasets from a pconf % pruned haplotype file
##############################################################################

library(data.table)
library(dplyr)

source("/home/math/siglersa/mastersProject/Input/create_haps_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/create_haps_funcs.R")


Pop = 'NFE'
# p_case = 100
# p_exp = 100
p_conf = 100
Nsim = 20000 
pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether
int_prune = 100
ext_prune = 80

# Number of individuals in each dataset
Ncase = Nint = 5000
Ncc = 10000 #cc10k = 10000, cc5k = 5000

# Haplotype column indices
cols = 1:40000

# Vectors to store the observed variants per MAC bin
# obs_MACbin_syn_pconf <- data.frame(matrix(ncol = 9, nrow = 100))
# colnames(obs_MACbin_syn_pconf) <- c('Singletons', 'Doubletons', 'MAC.3.5', 'MAC.6.10',
#                                   'MAC.11.20', 'MAC.21.MAF0.5', 'MAF0.5.1', 'rep', 'data')
# obs_MACbin_fun_pconf = obs_MACbin_fun_exp = obs_MACbin_syn_exp = obs_MACbin_syn_pconf
# obs_MACbin_fun_pcase = obs_MACbin_fun_exp = obs_MACbin_syn_pconf = obs_MACbin_syn_exp


mac_dir = '/home/math/siglersa/mastersProject/Input/'
# dir_in = '/storage/math/projects/compinfo/simulations/output/20K_NFE/'
# dir_in = '/storage/math/projects/RAREsim/Cases/Sim_20k/NFE/data/' #For pruning OG hap file
# dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/pruneDown/100v', p_conf, '/')
dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', int_prune, 'v', ext_prune, '/') #For pruning OG hap file
# dir_out = paste0('/home/math/siglersa/mastersProject/20K_NFE/cc10k/100v', p_conf, '/')
# dir_out = paste0('/home/math/siglersa/mastersProject/20K_NFE/pruneDown/100v', p_conf, '/')
dir_out = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', int_prune, 'v', ext_prune, '/datasets/') #For pruning OG hap file

# mac_dir = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/input/'
# dir_in = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/input/'
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/output/'


### read in the expected number of functional and synonymous variants from RAREsim
# exp_fun_case = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_fun_', p_case,  '.txt'), header=T, sep='\t')
# exp_syn = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_syn_', p_exp,  '.txt'), header=T, sep='\t')
# exp_fun = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_fun_', p_exp,  '.txt'), header=T, sep='\t')
# exp_fun_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_fun_', p_conf, '.txt'), header=T, sep='\t')
# exp_syn_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, '_', Pop, '_syn_', p_conf, '.txt'), header=T, sep='\t')


set.seed(1) # Will be different for each replicate but same for each run
# j = 1
for(j in 1:100){
  
  # read in the legend file 
  # leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.legend'), header=T, sep='\t')
  # leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.copy.legend'), header=T, sep='\t') #For pruning OG hap file
  # leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.s_only.legend'), header=T, sep='\t') # prune sequentially
  # leg$row = 1:nrow(leg)
  
  # read in the haplotype file
  # hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.', p_case, 'syn.haps.gz'))
  # hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.controls.haps.gz')) #For pruning OG hap file
  # hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_conf, '.haps.gz')) #prune together
  hap = fread(paste0(dir_in, 'chr19.block37.', Pop, '.sim', j, '.all.', p_conf, 'fun.', p_conf, 'syn.haps.gz')) # prune separately, sequentially
  hap = as.data.frame(hap)
  
  # add allele counts to the haplotypes
  # leg$count = rowSums(hap)
  
  # convert to minor allele counts
  # leg$MAC = ifelse(leg$count>Nsim, 2*Nsim-leg$count, leg$count)
  
  # subset the legend file to the functional variants (those are the only ones we'll prune)
  # leg_fun = leg %>% filter(fun=="fun")
  # leg_syn = leg %>% filter(fun=="syn")
  
  # Select the columns for necessary data sets
  cases = sort(sample(x=cols, size = 2*Ncase, replace = FALSE))
  int = sort(sample(x=cols[! cols %in% cases], size = 2*Nint, replace = FALSE))
  cc = sort(sample(x=cols[! cols %in% c(cases, int)], size = 2*Ncc, replace = FALSE))
  
  # subset the case haplotypes (120% fun and 100% syn)
  # hap_cases_pcase = hap[, cases]
  
  # write the haplotype file for the cases (power)
  # fwrite(hap_cases_pcase, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.cases.', p_case, 'fun.100syn.haps.gz'),
  #        quote=F, row.names=F, col.names=F, sep=' ')
  
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
  
  # PRUNE OG Hap file to 100% fun and 100% syn
  # rem_fun = select_var(leg_fun, exp_fun)
  # rem_syn = select_var(leg_syn, exp_syn)
  # hap_exp = prune_var(rbind(rem_fun, rem_syn), hap, Nsim)
  
  # subset the pruned haplotypes for 100% pruned
  # hap_cases = hap_all_pruned[, cases]
  # hap_int = hap_all_pruned[, int]
  # hap_cc = hap_all_pruned[, cc]
  # hap_cases = hap_exp[, cases]
  # hap_int = hap_exp[, int]
  # hap_cc = hap_exp[, cc]
  hap_cases = hap[, cases]
  hap_int = hap[, int]
  hap_cc = hap[, cc]
  
  # hap_all_pruned = hap # just need if going from 100% to 80%
  
  # write the haplotype files for the cases (type I error), internal and common controls
  fwrite(hap_cases, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')

  fwrite(hap_int, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')

  fwrite(hap_cc, paste0(dir_out, 'chr19.block37.', Pop, '.sim', j, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
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
  
  print(j)
}

# Save the obs MAC bin files
# fwrite(obs_MACbin_fun_pcase, paste0(dir_out, 'obs_MACbin_fun', p_case, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_fun_exp, paste0(dir_out, 'obs_MACbin_fun', p_exp, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_fun_pconf, paste0(dir_out, 'obs_MACbin_fun', p_conf, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_syn_exp, paste0(dir_out, 'obs_MACbin_syn', p_exp, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# fwrite(obs_MACbin_syn_pconf, paste0(dir_out, 'obs_MACbin_syn', p_conf, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
