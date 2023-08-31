############################################################################## 
# This file is used to test type I error for the three method with the additional
# scenario of using proxECAT and LogProx to test cases vs internal controls
# Mainly used for testing 100% vs 100% and 80% vs 80% pruned files for a
# homogeneous population
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

# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/read_in_funcs.R")
# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/general_data_manip.R")
# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/methods_funcs.R")
# source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/create_haps_funcs.R")

Pop = 'NFE'
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 99
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)


# Set appropriate directories
dir_leg = '/storage/math/projects/compinfo/simulations/output/20K_NFE/'
dir_in = '/home/math/siglersa/mastersProject/20K_NFE/cc10k/'
dir_out ='/home/math/siglersa/mastersProject/Results/cc10k/'

# dir_leg = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'
# dir_in = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'

# create empty vectors to store the p-values from each replicate
prox_p = prox_int_p = c()
prox2_p = prox2_all_p = prox2_int_p =  c()
iecat_p = c()
skat_int_p = skat_ext_p = skat_all_p = c()

# Create dataframe to store counts and ratios of fun:syn alleles for each dataset 
# ratios <- data.frame(matrix(ncol = 4, nrow = 3))
# colnames(ratios) <- c('Dataset', 'Functional', 'Synonymous', 'Ratio')
# ratios[, "Dataset"] <- c("Cases", "Internal Controls", "External Controls")


set.seed(1) 
# i=1
# loop through the simulation replicates
for (i in 1:100){
  
  # read in the legend file
  leg = read_leg_homo(dir_leg, Pop, i)
  
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
  
  ### CHECK COUNTS BEFORE PROXECAT
  # counts.prox = c()
  # 
  # case.fun = rare_counts(count_cases, leg_fun, leg_syn, maf)
  # counts.prox = c(counts.prox, c(case.fun[1], case.fun[2]))
  # 
  # ctrl.fun = rare_counts(count_int, leg_fun, leg_syn, maf)
  # counts.prox = c(counts.prox, c(ctrl.fun[1], ctrl.fun[2]))
  # 
  # prox = proxecat(counts.prox[1], counts.prox[2], counts.prox[3], counts.prox[4])
  ################################
  
  ### Run proxECAT and extract p-value 
  prox = prox_data_prep(leg_fun, leg_syn, count_cases, count_cc, maf)
  prox_int = prox_data_prep(leg_fun, leg_syn, count_cases, count_int, maf)
  
  # store proxECAT p-values
  prox_p = c(prox_p, prox)
  prox_int_p = c(prox_int_p, prox_int)
  
  ### Run LogProx and extract p-value
  p_prox2 = logprox_data_prep(leg, count_cases, count_int, count_cc, common_ext, common_all, adj=FALSE)
  p_prox2_int = logprox_int_prep(leg, count_cases, count_int, common_int)
  
  # NEED TO CHECK INDEXING
  # Store LogProx p-values
  prox2_p = c(prox2_p, p_prox2[[1]])
  prox2_all_p = c(prox2_all_p, p_prox2[[2]])
  prox2_int_p = c(prox2_int_p, p_prox2_int)
  
  ### Run iECAT and extract p-value
  run_iecat = iecat_data_prep(geno_cases, geno_int, leg, common_all, count_cc, Ncc)
  
  # Store iECAT p-values
  iecat_p = c(iecat_p, run_iecat[[1]])
  
  # Run SKAT-O and extract p-values
  run_skat = skat_data_prep(geno_cases, geno_int, geno_cc, leg, common_ext, common_all)
  
  # Store SKAT p-values
  skat_int_p = c(skat_int_p, run_iecat[[2]])
  skat_ext_p = c(skat_ext_p, run_skat[[1]])
  skat_all_p = c(skat_all_p, run_skat[[2]])
  
  print(i)
}

# Combine p-values from each method into one dataframe
results = data.frame(prox_p, prox_int_p, prox2_p, prox2_all_p, prox2_int_p, 
                     iecat_p, skat_int_p, skat_ext_p, skat_all_p)

# Save results
write.table(results, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F)
