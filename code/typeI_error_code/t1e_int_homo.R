############################################################################## 
# This file is used to test type I error for the three method with the additional
# scenario of using proxECAT and LogProx to test cases vs internal controls
# on a homogeneous population
##############################################################################

# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(devtools)
library(proxecat)
library(SKAT)
library(iECAT)
library(Summix)

source("/home/math/siglersa/mastersProject/Input/read_in_funcs.R")
source("/home/math/siglersa/mastersProject/Input/general_data_manip.R")
source("/home/math/siglersa/mastersProject/Input/methods_funcs.R")
source("/home/math/siglersa/mastersProject/Input/create_haps_funcs.R")

# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/pruning_code/create_haps_funcs.R")

Pop = 'NFE'
pruning = 'pruneSequentially' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim
folder = '100v100' # For testing 100v100, 100v80, and 80v80 from the 100v80 data
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 100
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)


# Set appropriate directories
# dir_leg = '/storage/math/projects/compinfo/simulations/output/20K_NFE/'
# dir_leg = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', int_prune, 'v', ext_prune, '/')
# dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', int_prune, 'v', ext_prune, '/')
dir_leg = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/')
dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/')
# dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/datasets/')
# dir_out ='/home/math/siglersa/mastersProject/Results/cc10k/'
dir_out = paste0('/home/math/siglersa/mastersProject/Output/', pruning, '/', folder, '/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_leg = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/input/pruneTogether/'
# dir_in = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/input/pruneTogether/'
# dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/'

# create empty vectors to store the p-values from each replicate
prox_p = prox_int_p = c()
prox2_p = prox2_all_p = prox2_int_p =  c()
iecat_p = c()
skat_int_p = skat_ext_p = skat_all_p = c()
# skat_int_p_fun = skat_ext_p_fun = skat_all_p_fun = c()
iecat_p_syn = c()
skat_int_p_syn = skat_ext_p_syn = skat_all_p_syn = c()

# Create dataframe to store counts and ratios of fun:syn alleles for each dataset 
# ratios <- data.frame(matrix(ncol = 4, nrow = 13))
# colnames(ratios) <- c('Dataset', 'Functional', 'Synonymous', 'Ratio')
# ratios[, "Dataset"] <- c("Cases-Hap", "Internal Controls-Hap", "External Controls-Hap",
#                          "Cases-ProxECAT Internal", "Internal Controls-ProxECAT",
#                          "Cases-ProxECAT External", "External Controls-ProxECAT",
#                          "Cases-LogProx Internal", "Internal Controls-LogProx Internal",
#                          "Cases-LogProx External", "External Controls-LogProx External",
#                          "Cases-Logprox Internal + External", "Controls-Logprox Internal + External")
# 
# pvals <- data.frame((matrix(ncol = 2, nrow = 5)))
# colnames(pvals) <- c('Method', 'P-value')
# pvals[, "Method"] <- c('ProxECAT-Internal', 'ProxECAT-External', 'LogProx-Internal',
#                        'LogProx-External', 'LogProx-Internal+External')

# proxEcounts <- data.frame(matrix(ncol = 11, nrow = 100))
# colnames(proxEcounts) <- c('Case-Fun (O)', 'Case-Syn (O)', 'Control-Fun (O)', 'Control-Syn (O)', 
#                            'Control-Fun (E)', 'Control-Syn (E)', 'Control-Fun (O-E)',
#                            'Control-Syn (O-E)', 'Ratio-Case', 'Ratio-Control', 'P-Value')


set.seed(1) 
# i=1
# loop through the simulation replicates
for (i in 1:100){
  
  # read in the legend file
  leg = read_leg_homo(dir_leg, Pop, i)
  # leg = read.table(paste0(dir_in, 'chr19.block37.', Pop, '.sim', i, '.s_only.legend'), header=T, sep='\t') # prune sequentially
  # leg$row = 1:nrow(leg) # prune sequentially
  
  leg_fun = leg %>% filter(fun=="fun")
  leg_syn = leg %>% filter(fun=="syn")
  
  # read in the haplotype and reference files
  hap_cases = read_hap_homo(dir_in, Pop, i, "cases", p_case_fun, p_case_syn)
  hap_int = read_hap_homo(dir_in, Pop, i, "internal.controls", p_int_fun, p_int_syn)
  hap_cc = read_hap_homo(dir_in, Pop, i, "common.controls", p_cc_fun, p_cc_syn)
  
  ### CHECK COUNTS
  ### Check ratio of fun to syn rare alleles in cases
  # ratios[1, 2] = rare_var(leg_fun$row, hap_cases, maf = maf)
  # ratios[1, 3] = rare_var(leg_syn$row, hap_cases, maf = maf)
  # ratios[1, 4] = ratios[1, 2]/ratios[1, 3]
  # 
  # ## Check ratio of fun to syn variants in internal controls
  # ratios[2, 2] = rare_var(leg_fun$row, hap_int, maf = maf)
  # ratios[2, 3] = rare_var(leg_syn$row, hap_int, maf = maf)
  # ratios[2, 4] = ratios[2, 2]/ratios[2, 3]
  # 
  # ### Check ratio of fun to syn variants in external controls
  # ratios[3, 2] = rare_var(leg_fun$row, hap_cc, maf = maf)
  # ratios[3, 3] = rare_var(leg_syn$row, hap_cc, maf = maf)
  # ratios[3, 4] = ratios[3, 2]/ratios[3, 3]
  ################
  
  # convert the haplotypes into genotypes
  geno_cases = make_geno(hap_cases)
  geno_int = make_geno(hap_int)
  geno_cc = make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_cases = calc_allele_freqs(geno_cases, Ncase)
  count_int = calc_allele_freqs(geno_int, Nint)
  count_cc = calc_allele_freqs(geno_cc, Ncc)
  count_all = calc_allele_freqs_all(count_cases, count_int, count_cc, Ncase, Nint, Ncc)
  
  # identify the common variants
  common_ext = leg[which(count_cases$maf > maf | count_cc$maf > maf),]
  common_all = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc$maf > maf),]
  common_int = leg[which(count_cases$maf > maf | count_int$maf > maf),]
  
  ### Run proxECAT and extract p-value 
  # prox = prox_data_prep(leg_fun, leg_syn, count_cases, count_cc, maf)
  # prox_int = prox_data_prep(leg_fun, leg_syn, count_cases, count_int, maf)
  prox = prox_data_prep(leg, count_cases, count_cc, common_ext, adj=FALSE)
  prox_int = prox_int_prep(leg, count_cases, count_int, common_int)
  
  # store proxECAT p-values
  prox_p = c(prox_p, prox)
  prox_int_p = c(prox_int_p, prox_int)

  ##############################################################################
  ### Checks
  # convert genotypes into long format for ProxECAT v2
  # data.cases = make_long(count_cases, leg, "case", "int")
  # data.cc = make_long(count_cc, leg, "control", "ext")
  # data.int = make_long(count_int, leg, "control", "int")
  # 
  # # combine the data together AND REMOVE COMMON VARIANTS
  # data.prox = data.frame(lapply(rbind(data.cases, data.cc), factor)) %>%
  #   filter(!(id %in% common_ext$id))
  # 
  # data.prox.int = data.frame(lapply(rbind(data.cases, data.int), factor)) %>%
  #   filter(!(id %in% common_int$id))
  # 
  # # getting overall counts for functional & case status
  # # data for proxECAT method
  # counts.prox = data.prox %>% count(case, fun)
  # counts.prox.int = data.prox.int %>% count(case, fun)
  # 
  # # Run proxECAT
  # prox = proxecat(counts.prox$n[1], counts.prox$n[2], counts.prox$n[3], counts.prox$n[4])
  # prox.int = proxecat(counts.prox.int$n[1], counts.prox.int$n[2], counts.prox.int$n[3], counts.prox.int$n[4])
  # 
  # pvals[1:2, 2] <- c(prox.int$p.value, prox$p.value)
  
  ### Check ratio of fun to syn rare alleles in cases-Internal
  # ratios[4, 2:3] = c(counts.prox.int$n[1], counts.prox.int$n[2])
  # ratios[4, 4] = ratios[4, 2]/ratios[4, 3]
  # 
  # ## Check ratio of fun to syn variants in internal controls
  # ratios[5, 2:3] = c(counts.prox.int$n[3], counts.prox.int$n[4])
  # ratios[5, 4] = ratios[5, 2]/ratios[5, 3]
  # 
  # ### Check ratio of fun to syn rare alleles in cases-External
  # ratios[6, 2:3] = c(counts.prox$n[1], counts.prox$n[2])
  # ratios[6, 4] = ratios[6, 2]/ratios[6, 3]
  # 
  # ### Check ratio of fun to syn variants in external controls
  # ratios[7, 2:3] = c(counts.prox$n[3], counts.prox$n[4])
  # ratios[7, 4] = ratios[7, 2]/ratios[7, 3]
  
  # For proxecat expanded
  # proxEcounts[i, 1:4] <- counts.prox # store counts
  # proxEcounts[i, 'Ratio-Case'] <- proxEcounts[i, 1]/proxEcounts[i, 2] # calc case ratios
  # proxEcounts[i, 'Ratio-Control'] <- proxEcounts[i, 3]/proxEcounts[i, 4] # calc ctrl ratios
  
  # For testing 100% v ext_prune
  # proxEcounts[i, 'Control-Fun (E)'] <- proxEcounts[i, 'Case-Fun (O)']*2*(ext_prune/100) # calc E ctrl-fun
  # proxEcounts[i, 'Control-Syn (E)'] <- proxEcounts[i, 'Case-Syn (O)']*2*(ext_prune/100) # calc E ctrl-syn
  # # For testing ext_prune v ext_prune
  # proxEcounts[i, 'Control-Fun (E)'] <- proxEcounts[i, 'Case-Fun (O)']*2 # calc E ctrl-fun
  # proxEcounts[i, 'Control-Syn (E)'] <- proxEcounts[i, 'Case-Syn (O)']*2 # calc E ctrl-syn
  # 
  # proxEcounts[i, 'Control-Fun (O-E)'] <- proxEcounts[i, 'Control-Fun (O)']-proxEcounts[i, 'Control-Fun (E)'] # calc O-E ctrl-fun
  # proxEcounts[i, 'Control-Syn (O-E)'] <- proxEcounts[i, 'Control-Syn (O)']-proxEcounts[i, 'Control-Syn (E)'] # calc O-E ctrl-syn
  
  # Run proxECAT
  # prox = proxecat(counts.prox[1], counts.prox[2], counts.prox[3], counts.prox[4])
  # proxEcounts[i, 'P-Value'] <- prox$p.value
  ##############################################################################
  
  ### Run LogProx and extract p-value
  p_prox2 = logprox_data_prep(leg, count_cases, count_int, count_cc, common_ext, common_all, adj=FALSE)
  p_prox2_int = logprox_int_prep(leg, count_cases, count_int, common_int)
  
  # Store LogProx p-values
  prox2_p = c(prox2_p, p_prox2[[1]])
  prox2_all_p = c(prox2_all_p, p_prox2[[2]])
  prox2_int_p = c(prox2_int_p, p_prox2_int)
  
  ##############################################################################
  #### Check LogProx Counts
  # convert genotypes into long format for ProxECAT v2
  # data.cases = make_long(count_cases, leg, "case", "int")
  # data.int = make_long(count_int, leg, "control", "int")
  # data.cc = make_long(count_cc, leg, "control", "ext")
  # 
  # # combine the data together AND REMOVE COMMON VARIANTS
  # data.prox = data.frame(lapply(rbind(data.cases, data.cc), factor)) %>% 
  #   filter(!(id %in% common_ext$id))
  # 
  # data.int = data.frame(lapply(rbind(data.cases, data.int), factor)) %>% 
  #   filter(!(id %in% common_int$id))
  # 
  # data.all = data.frame(lapply(rbind(data.cases, data.int, data.cc), factor)) %>% 
  #   filter(!(id %in% common_all$id))
  # 
  # counts.logprox = data.prox %>% count(case, fun)
  # counts.int.logprox = data.int %>% count(case, fun)
  # counts.all.logprox = data.all %>% count(case, fun)
  # 
  # glm.prox = glm(fun ~ case, data=data.prox, family="binomial") 
  # glm.int.prox = glm(fun ~ case, data=data.int, family="binomial")
  # glm.all.prox = glm(fun ~ case + group, data=data.all, family="binomial")
  
  # save the p-value for case/control status
  # p.prox = summary(glm.prox)$coefficients[2,4]
  # p.prox.int = summary(glm.int.prox)$coefficients[2,4]
  # p.prox.all = summary(glm.all.prox)$coefficients[2,4]
  # 
  # pvals[3:5, 2] <- c(p.prox.int, p.prox, p.prox.all)
  
  ### Check ratio of fun to syn rare alleles in cases-internal
  # ratios[8, 2:3] = counts.int.logprox$n[1:2]
  # ratios[8, 4] = ratios[8, 2]/ratios[8, 3]
  # 
  # ## Check ratio of fun to syn variants in internal controls
  # ratios[9, 2:3] = counts.int.logprox$n[3:4]
  # ratios[9, 4] = ratios[9, 2]/ratios[9, 3]
  # 
  # ### Check ratio of fun to syn variants in cases-external
  # ratios[10, 2:3] = counts.logprox$n[1:2]
  # ratios[10, 4] = ratios[10, 2]/ratios[10, 3]
  # 
  # ### Check ratio of fun to syn variants in external controls
  # ratios[11, 2:3] = counts.logprox$n[3:4]
  # ratios[11, 4] = ratios[11, 2]/ratios[11, 3]
  # 
  # ### Check ratio of fun to syn rare alleles in cases-internal+external
  # ratios[12, 2:3] = counts.all.logprox$n[1:2]
  # ratios[12, 4] = ratios[12, 2]/ratios[12, 3]
  # 
  # ## Check ratio of fun to syn variants in controls-internal+external
  # ratios[13, 2:3] = counts.all.logprox$n[3:4]
  # ratios[13, 4] = ratios[13, 2]/ratios[13, 3]
  # 
  # write.csv(ratios, paste0(dir_out, "fixed_proxECAT_v_logProx_checks_sim1_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
  # write.csv(pvals, paste0(dir_out, "fixed_proxECAT_v_logProx_pvalues_sim1_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
  ##############################################################################
  
  ### Run iECAT and extract p-value
  # run_iecat = iecat_data_prep(geno_cases, geno_int, leg, common_all, count_cc, Ncc)
  run_iecat = iecat_data_prep(geno_cases, geno_int, leg_syn, common_all, count_cc, Ncc)
  run_iecat_syn = iecat_data_prep(geno_cases, geno_int, leg_fun, common_all, count_cc, Ncc)
  
  # Store iECAT p-values
  iecat_p = c(iecat_p, run_iecat[[1]])
  # iecat_p_fun = c(iecat_p_fun, run_iecat_fun[[1]])
  iecat_p_syn = c(iecat_p_syn, run_iecat_syn[[1]])
  
  # Run SKAT-O and extract p-values
  # run_skat = skat_data_prep(geno_cases, geno_int, geno_cc, leg, common_ext, common_all)
  run_skat = skat_data_prep(geno_cases, geno_int, geno_cc, leg_syn, common_ext, common_all)
  run_skat_syn = skat_data_prep(geno_cases, geno_int, geno_cc, leg_fun, common_ext, common_all)
  
  # Store SKAT p-values
  skat_int_p = c(skat_int_p, run_iecat[[2]])
  skat_ext_p = c(skat_ext_p, run_skat[[1]])
  skat_all_p = c(skat_all_p, run_skat[[2]])
  # skat_int_p_fun = c(skat_int_p_fun, run_iecat_fun[[2]])
  # skat_ext_p_fun = c(skat_ext_p_fun, run_skat_fun[[1]])
  # skat_all_p_fun = c(skat_all_p_fun, run_skat_fun[[2]])
  skat_int_p_syn = c(skat_int_p_syn, run_iecat_syn[[2]])
  skat_ext_p_syn = c(skat_ext_p_syn, run_skat_syn[[1]])
  skat_all_p_syn = c(skat_all_p_syn, run_skat_syn[[2]])
  
  print(i)
}

# Combine p-values from each method into one dataframe
# results = data.frame(prox_p, prox_int_p, prox2_p, prox2_all_p, prox2_int_p,
#                      iecat_p, skat_int_p, skat_ext_p, skat_all_p)

# results = data.frame(iecat_p_syn, skat_int_p_syn, skat_ext_p_syn, skat_all_p_syn)

results = data.frame(prox_p, prox_int_p, prox2_p, prox2_all_p, prox2_int_p,
                     iecat_p, skat_int_p, skat_ext_p, skat_all_p,
                     iecat_p_syn, skat_int_p_syn, skat_ext_p_syn, skat_all_p_syn)

# Save results
# write.table(results, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F)
write.table(results, paste0(dir_out, "T1e_",  pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F)
# fwrite(proxEcounts, paste0(dir_out, 'proxECAT_counts_expanded_', Pop, '_', int_prune, "_v_", ext_prune, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
