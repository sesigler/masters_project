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


Pop1 = 'AFR'
Pop2 = 'NFE'
scen = 's1' #scenario = 's1' or 's2'
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 100
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
Nref = 500
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
pi_tar1 = 0.15 #pi.target for NFE: 0.15 for s1 or 0 for s2
pi_tar2 = 0.85 #pi.target for AFR: 0.85 for s1 or 1 for s2


# Set appropriate directories
dir_leg ='/storage/math/projects/compinfo/simulations/output/NFE_AFR_pops/'
dir_in = '/home/math/siglersa/mastersProject/new_AFR_NFE_pops/cc10k/'
dir_out ='/home/math/siglersa/mastersProject/Results/cc10k/'


# create empty vectors to store the p-values from each replicate
prox_p = prox_p_adj = c()
prox2_p = prox2_all_p = prox2_p_adj = prox2_all_p_adj = c()
iecat_p = iecat_p_adj = c()
## OPTIONAL
prox_int_p = prox_int_p_adj = prox2_int_p = prox2_int_p_adj = c() 
# Store summix ancestry proportion estimates
prop_ests_cc = prop_ests_cases = prop_ests_int = c()


set.seed(1) 
# loop through the simulation replicates
for (i in 1:100){
  
  # read in the legend file
  leg = read_leg(dir_leg, Pop1, Pop2, i)

  # read in the haplotype and reference files
  hap_cases = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn)
  hap_int = read_hap(dir_in, Pop1, Pop2, i, scen, "internal.controls", p_int_fun, p_int_syn)
  hap_cc = read_hap(dir_in, Pop1, Pop2, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  hap_ref_afr = read_ref(dir_in, Pop1, i, scen)
  hap_ref_nfe = read_ref(dir_in, Pop2, i, scen)

  # convert the haplotypes into genotypes
  geno_cases = make_geno(hap_cases)
  geno_int = make_geno(hap_int)
  geno_cc = make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_cases = calc_allele_freqs(geno_cases, Ncase)
  count_int = calc_allele_freqs(geno_int, Nint)
  count_cc = calc_allele_freqs(geno_cc, Ncc)
  
  count_all = calc_allele_freqs_all(count_cases, count_int, count_cc, Ncase, Nint, Ncc)
  count_ref_afr = calc_allele_freqs_ref(Pop1, hap_ref_afr, Nref)
  count_ref_nfe = calc_allele_freqs_ref(Pop2, hap_ref_nfe, Nref)
  
  cc_refs = cbind(count_cc, count_ref_afr, count_ref_nfe)
  cases_refs = cbind(count_cases, count_ref_afr, count_ref_nfe)
  int_refs = cbind(count_int, count_ref_afr, count_ref_nfe)

  # Estimate ancestry proportions using only COMMON variants
  cc_est_prop = est_props(cc_refs, maf)
  cases_est_prop = est_props(cases_refs, maf)
  int_est_prop = est_props(int_refs, maf)
  
  prop_ests_cc <- rbind(prop_ests_cc, cc_est_prop)
  prop_ests_cases <- rbind(prop_ests_cases, cases_est_prop)
  prop_ests_int <- rbind(prop_ests_int, int_est_prop)
  
  # Calculate adjusted AFs 
  count_cc_adj = calc_adjusted_AF(cc_refs, cc_est_prop, pi_tar1, pi_tar2, Ncc)
  
  # identify the common variants
  common_ext = leg[which(count_cases$maf > maf | count_cc$maf > maf),]
  common_all = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc$maf > maf),]
  
  common_ext_adj = leg[which(count_cases$maf > maf | count_cc_adj$maf > maf),]
  common_all_adj = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc_adj$maf > maf),]
  
  ### Run proxECAT and extract p-value 
  prox = prox_data_prep(leg, count_cases, count_cc, common_ext, adj=FALSE)
  prox_adj = prox_data_prep(leg, count_cases, count_cc_adj, common_ext_adj, adj=TRUE)
  
  # store proxECAT p-values
  prox_p = c(prox_p, prox)
  prox_p_adj = c(prox_p_adj, prox_adj)
  
  ### Run LogProx and extract p-value
  p_prox2 = logprox_data_prep(leg, count_cases, count_int, count_cc, common_ext, common_all, adj=FALSE)
  p_prox2_adj = logprox_data_prep(leg, count_cases, count_int, count_cc_adj, 
                                  common_ext_adj, common_all_adj, adj=TRUE)

  # NEED TO CHECK INDEXING
  # Store LogProx p-values
  prox2_p = c(prox2_p, p_prox2[[1]])
  prox2_all_p = c(prox2_all_p, p_prox2[[2]])
  
  prox2_p_adj = c(prox2_p_adj, p_prox2_adj[[1]])
  prox2_all_p_adj = c(prox2_all_p_adj, p_prox2_adj[[2]])
  
  ### Run iECAT and extract p-value
  run_iecat = iecat_data_prep(geno_cases, geno_int, leg, common_all, count_cc, Ncc)
  run_iecat_adj = iecat_data_prep(geno_cases, geno_int, leg, common_all_adj, count_cc_adj, Ncc)
  
  # Store iECAT p-values
  iecat_p = c(iecat_p, run_iecat)
  iecat_p_adj = c(iecat_p_adj, run_iecat_adj)
  
  print(i)
}

# Combine p-values from each method into one dataframe
results = data.frame(prox_p, prox2_p, prox2_all_p, iecat_p)
results_adj = data.frame(prox_p_adj, prox2_p_adj, prox2_all_p_adj, iecat_p_adj)

# Save results
write.table(results, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
write.table(results_adj, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)

# Proportion estimates
write.table(prop_ests_cases, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_case_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
write.table(prop_ests_cc, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_cc_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
write.table(prop_ests_int, paste0(dir_out, "T1e_", int_prune, "_v_", ext_prune, "_int_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), quote=F, row.names=F)
