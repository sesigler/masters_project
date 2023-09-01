############################################################################## 
# This file is used to test type I error for the three method with the additional
# scenario of using proxECAT and LogProx to test cases vs internal controls
# Mainly used for testing 100% vs 100% and 80% vs 80% pruned files
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

# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/read_in_funcs.R")
# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/general_data_manip.R")
# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/methods_funcs.R")
# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/create_haps_funcs.R")


Pop1 = 'AFR'
Pop2 = 'NFE'
scen = 's1' #scenario = 's1' or 's2'
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 99
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
Nref = 500
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
pi_tar1 = 0.15 #pi.target for NFE: 0.15 for s1 or 0 for s2
pi_tar2 = 0.85 #pi.target for AFR: 0.85 for s1 or 1 for s2


# Set appropriate directories
dir_leg ='/storage/math/projects/compinfo/simulations/output/NFE_AFR_pops/'
dir_in = '/home/math/siglersa/mastersProject/new_AFR_NFE_pops/cc10k/100v99/'
# dir_out ='/home/math/siglersa/mastersProject/Results/cc10k/'

# dir_in = '/home/math/siglersa/mastersProject/Input/'
dir_out = '/home/math/siglersa/mastersProject/Output/'


# dir_leg = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'
# dir_leg = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/AFR_NFE_pops/cc10K/100v99/'
# dir_in = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'
# dir_in = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/AFR_NFE_pops/cc10K/100v99/'
# dir_out = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/output/'

# create empty vectors to store the p-values from each replicate
# prox_p = prox_p_adj = c()
# prox2_p = prox2_all_p = prox2_p_adj = prox2_all_p_adj = c()
# iecat_p = iecat_p_adj = c()
# skat_int_p = skat_ext_p = skat_all_p = c()
# # Added for cases vs internal controls testing
# prox_int_p = prox2_int_p = c()

prox_p = prox2_p = c()

# Create dataframe to store counts and ratios of fun:syn alleles for each dataset 
# ratios <- data.frame(matrix(ncol = 4, nrow = 5))
# colnames(ratios) <- c('Dataset', 'Functional', 'Synonymous', 'Ratio')
# ratios[, "Dataset"] <- c("Cases", "Internal Controls", "External Controls", "Ref AFR", "Ref NFE")

proxEcounts <- data.frame(matrix(ncol = 6, nrow = 11))
colnames(proxEcounts) <- c('Case-Fun', 'Case-Syn', 'Control-Fun', 'Control-Syn', 'Ratio-Case', 'Ratio-Control')

logprox_counts <- data.frame(matrix(ncol = 6, nrow = 11))
colnames(logprox_counts) <- c('Case-Fun', 'Case-Syn', 'Control-Fun', 'Control-Syn', 'Ratio-Case', 'Ratio-Control')

# Read in csv file with rows to check
sim_rows <- read.csv(paste0(dir_out, 'rows_99_v_80_counts_to_check.csv'), header=T)
sim_reps <- unique(sim_rows$row)
i = sim_reps[1]
j = 1
set.seed(1) 
# loop through the simulation replicates
# for (i in 1:100){
for (i in sim_reps){
  
  # read in the legend file
  leg = read_leg(dir_leg, Pop1, Pop2, i)
  
  leg_fun = leg %>% filter(fun=="fun")
  leg_syn = leg %>% filter(fun=="syn")
  
  # read in the haplotype and reference files
  hap_cases = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn)
  hap_int = read_hap(dir_in, Pop1, Pop2, i, scen, "internal.controls", p_int_fun, p_int_syn)
  hap_cc = read_hap(dir_in, Pop1, Pop2, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  hap_ref_afr = read_ref(dir_in, Pop1, i, scen)
  hap_ref_nfe = read_ref(dir_in, Pop2, i, scen)
  
  ### CHECK COUNTS AND OUTPUT TO CSV
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
  # 
  # ### Check ratio of fun to syn variants in ref AFR
  # ratios[4, 2] = rare_var(leg_fun$row, hap_ref_afr, maf = maf)
  # ratios[4, 3] = rare_var(leg_syn$row, hap_ref_afr, maf = maf)
  # ratios[4, 4] = ratios[4, 2]/ratios[4, 3]
  # 
  # ### Check ratio of fun to syn variants in ref NFE
  # ratios[5, 2] = rare_var(leg_fun$row, hap_ref_nfe, maf = maf)
  # ratios[5, 3] = rare_var(leg_syn$row, hap_ref_nfe, maf = maf)
  # ratios[5, 4] = ratios[5, 2]/ratios[5, 3]
  # 
  # # Output ratios
  # fwrite(ratios, paste0(dir_out, 'ratios_100_v_',  ext_prune, '_sim', i, '_t1e_checks.csv'),
  #        quote=F, row.names=F, col.names=T, sep=',')
  ##################################
  
  # convert the haplotypes into genotypes
  geno_cases = make_geno(hap_cases)
  geno_int = make_geno(hap_int)
  geno_cc = make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_cases = calc_allele_freqs(geno_cases, Ncase)
  count_int = calc_allele_freqs(geno_int, Nint)
  count_cc = calc_allele_freqs(geno_cc, Ncc)
  
  # count_all = calc_allele_freqs_all(count_cases, count_int, count_cc, Ncase, Nint, Ncc)
  # count_ref_afr = calc_allele_freqs_ref(Pop1, hap_ref_afr, Nref)
  # count_ref_nfe = calc_allele_freqs_ref(Pop2, hap_ref_nfe, Nref)
  # 
  # cc_refs = cbind(count_cc, count_ref_afr, count_ref_nfe)
  
  # Estimate ancestry proportions using only COMMON variants
  # cc_est_prop = est_props(cc_refs, Pop1, Pop2, maf)
  
  # Calculate adjusted AFs 
  # count_cc_adj = calc_adjusted_AF(cc_refs, Pop2, cc_est_prop, pi_tar1, pi_tar2, Ncc)
  
  # identify the common variants
  common_ext = leg[which(count_cases$maf > maf | count_cc$maf > maf),]
  # common_all = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc$maf > maf),]
  
  # common_ext_adj = leg[which(count_cases$maf > maf | count_cc_adj$maf > maf),]
  # common_all_adj = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc_adj$maf > maf),]
  # 
  # common_int = leg[which(count_cases$maf > maf | count_int$maf > maf),]
  
  ### Run proxECAT and extract p-value 
  # prox = prox_data_prep(leg_fun, leg_syn, count_cases, count_cc, maf)
  counts.prox = c()
  
  case.fun = rare_counts(count_cases, leg_fun, leg_syn, maf)
  counts.prox = c(counts.prox, c(case.fun[1], case.fun[2]))
  
  ctrl.fun = rare_counts(count_cc, leg_fun, leg_syn, maf)
  counts.prox = c(counts.prox, c(ctrl.fun[1], ctrl.fun[2]))
  
  proxEcounts[j, 1:4] <- counts.prox
  proxEcounts[j, 5] <- proxEcounts[j, 1]/proxEcounts[j, 2]
  proxEcounts[j, 6] <- proxEcounts[j, 3]/proxEcounts[j, 4]
  
  # Run proxECAT
  prox = proxecat(counts.prox[1], counts.prox[2], counts.prox[3], counts.prox[4])
  # prox_adj = prox_data_prep(leg_fun, leg_syn, count_cases, count_cc_adj, maf)
  # prox_int = prox_data_prep(leg_fun, leg_syn, count_cases, count_int, maf)
  
  # store proxECAT p-values
  prox_p = c(prox_p, prox$p.value)
  # prox_p_adj = c(prox_p_adj, prox_adj)
  # prox_int_p = c(prox_int_p, prox_int)
  
  ### Run LogProx and extract p-value
  # p_prox2 = logprox_data_prep(leg, count_cases, count_int, count_cc, common_ext, common_all, adj=FALSE)
  
  data.cases = make_long(count_cases, leg, "case", "int")
  data.cc = make_long(count_cc, leg, "control", "ext")
  
  # combine the data together AND REMOVE COMMON VARIANTS
  data.prox = data.frame(lapply(rbind(data.cases, data.cc), factor)) %>% 
    filter(!(id %in% common_ext$id))
  
  logprox_counts[j, 1] <- length(which(data.prox$fun == "fun" & data.prox$case == "case"))
  logprox_counts[j, 2] <- length(which(data.prox$fun == "syn" & data.prox$case == "case"))
  logprox_counts[j, 3] <- length(which(data.prox$fun == "fun" & data.prox$case == "control"))
  logprox_counts[j, 4] <- length(which(data.prox$fun == "syn" & data.prox$case == "control"))
  
  logprox_counts[j, 5] <- logprox_counts[j, 1]/logprox_counts[j, 2]
  logprox_counts[j, 6] <- logprox_counts[j, 3]/logprox_counts[j, 4]
  # fit the ProxECATv2 model
  glm.prox = glm(fun ~ case, data=data.prox, family="binomial") 
  
  # save the p-value for case/control status
  p.prox = summary(glm.prox)$coefficients[2,4]
  # p.prox.all = summary(glm.all.prox)$coefficients[2,4]
  # p_prox2_adj = logprox_data_prep(leg, count_cases, count_int, count_cc_adj, 
  #                                 common_ext_adj, common_all_adj, adj=TRUE)
  # p_prox2_int = logprox_int_prep(leg, count_cases, count_int, common_int)
  
  # Store LogProx p-values
  prox2_p = c(prox2_p, p.prox)
  j <- j + 1
  # prox2_p = c(prox2_p, p_prox2[[1]])
  # prox2_all_p = c(prox2_all_p, p_prox2[[2]])
  # 
  # prox2_p_adj = c(prox2_p_adj, p_prox2_adj[[1]])
  # prox2_all_p_adj = c(prox2_all_p_adj, p_prox2_adj[[2]])
  # 
  # prox2_int_p = c(prox2_int_p, p_prox2_int)
  
  ### Run iECAT and extract p-value
  # run_iecat = iecat_data_prep(geno_cases, geno_int, leg, common_all, count_cc, Ncc)
  # run_iecat_adj = iecat_data_prep(geno_cases, geno_int, leg, common_all_adj, count_cc_adj, Ncc)
  # 
  # # Store iECAT p-values
  # iecat_p = c(iecat_p, run_iecat[[1]])
  # iecat_p_adj = c(iecat_p_adj, run_iecat_adj[[1]])
  # 
  # # Run SKAT-O and extract p-values
  # run_skat = skat_data_prep(geno_cases, geno_int, geno_cc, leg, common_ext, common_all)
  # 
  # # Store SKAT p-values
  # skat_int_p = c(skat_int_p, run_iecat[[2]])
  # skat_ext_p = c(skat_ext_p, run_skat[[1]])
  # skat_all_p = c(skat_all_p, run_skat[[2]])
  print(j)
  print(i)
}

# Combine p-values from each method into one dataframe
results = data.frame(prox_p, prox2_p)
# results = data.frame(prox_p, prox_int_p, prox2_p, prox2_all_p, prox2_int_p, 
#                      iecat_p, skat_int_p, skat_ext_p, skat_all_p)
# results_adj = data.frame(prox_p_adj, prox2_p_adj, prox2_all_p_adj, iecat_p_adj)

# Save results
write.table(results, paste0(dir_out, "T1e_checks_", int_prune, "_v_", ext_prune, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
fwrite(proxEcounts, paste0(dir_out, 'proxECAT_counts_', int_prune, "_v_", ext_prune, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
fwrite(logprox_counts, paste0(dir_out, 'logProx_counts_', int_prune, "_v_", ext_prune, '.csv'), quote=F, row.names=F, col.names=T, sep=',')
# write.table(results, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
# write.table(results_adj, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
