############################################################################## 
# This file is used to generate the haplotype files necessary for running
# the type I error and power calculations for proxECAT, LogProx, and iECAT-O
# on an admixed population
##############################################################################

library(data.table)
library(dplyr)

source("/home/math/siglersa/mastersProject/Input/create_haps_funcs.R")
source("/home/math/siglersa/mastersProject/Input/subset_haps_funcs.R")

source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/create_haps_funcs.R")
source("C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/subset_haps_funcs.R")


Pop1 = 'AFR'
Pop2 = 'NFE'
p_case = 120
p_conf = 99
Nsim = 22500 
maf = 0.001
scen = 's1' #scenario: 's1' or 's2'
Ncc = "cc10k"

# The following have all been multiplied by 2 for number of haplotypes
afr_case_size = afr_ic_size = 8500 #s1 = 8500, s2 = 10000 
nfe_case_size = nfe_ic_size = 1500 #s1 = 1500, s2 = 0 

afr_cc_size = 17000 #10kcc = 17000, 5kcc = 8500
nfe_cc_size = 3000 #10kcc = 3000, 5kcc = 1500

afr_ref_size = nfe_ref_size = 1000

# Haplotype column indices for each pop
afr_cols = 1:38000 
nfe_cols = 38001:45000

### Create empty vectors to store counts and ratios of fun:syn variants for each dataset
ratios_case_pow = ratios_case_t1e = ratios_case_pconf = ratios_int = ratios_int_pconf = ratios_cc = ratios_cc_pconf = c()


mac_dir = '/home/math/siglersa/mastersProject/Input/'
dir_in = '/storage/math/projects/compinfo/simulations/output/NFE_AFR_pops/'
dir_out = '/home/math/siglersa/mastersProject/new_AFR_NFE_pops/cc10k/'

mac_dir = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'
dir_in = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/'
dir_out = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/'


### read in the expected number of functional and synonymous variants from RAREsim
exp_fun_case = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_fun_', p_case,  '.txt'), header=T, sep='\t')
exp_syn = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_syn_100.txt'), header=T, sep='\t')
exp_fun = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_fun_100.txt'), header=T, sep='\t')
exp_fun_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_fun_', p_conf, '.txt'), header=T, sep='\t')
exp_syn_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_syn_', p_conf, '.txt'), header=T, sep='\t')


set.seed(1) # Will be different for each replicate but same for each run
j=1
for(j in 1:100){
  
  # read in the legend file
  leg = read.table(paste0(dir_in, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.legend'), header=T, sep='\t')
  leg$row = 1:nrow(leg)
  
  # read in the haplotype file
  hap = fread(paste0(dir_in, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.all.', p_case, 'fun.100syn.haps.gz'))
  hap = as.data.frame(hap)
  
  # add allele counts to the haplotypes
  leg$count = rowSums(hap)
  
  # convert to minor allele counts
  leg$MAC = ifelse(leg$count>Nsim, 2*Nsim-leg$count, leg$count)
  
  # subset the legend file to the functional variants (those are the only ones we'll prune)
  leg_fun = leg %>% filter(fun=="fun")
  leg_syn = leg %>% filter(fun=="syn")
  
  # Select AFR and NFE columns for necessary data sets
  cases = case_cols(afr_cols, nfe_cols, afr_case_size, nfe_case_size, scen)
  ics = int_cols(afr_cols, nfe_cols, afr_ic_size, nfe_ic_size, cases, scen)
  ccs = cc_cols(afr_cols, nfe_cols, afr_cc_size, nfe_cc_size, cases, ics, scen)
  refs = ref_cols(afr_cols, nfe_cols, afr_ref_size, nfe_ref_size, cases, ics, ccs, scen)
  
  # Create dataframe to store counts and ratios of fun:syn variants for each dataset 
  ratios <- data.frame(matrix(ncol = 4, nrow = 7))
  colnames(ratios) <- c('Dataset', 'Functional', 'Synonymous', 'Ratio')
  ratios[, "Dataset"] <- c("Cases (Power)", "Cases (T1E)", "Internal Controls", 
                           "External Controls", "External Controls (Pruned)",
                           "Cases (Pruned)", "Internal Controls (Pruned)")
  
  # subset the case haplotypes (120% fun and 100% syn)
  hap_cases_pcase = sub_hap_scen(hap, cases, scen)
  
  # write the haplotype file for the cases (power)
  fwrite(hap_cases_pcase, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.cases.', p_case, 'fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  ### Check ratio of fun to syn variants in cases (power)
  ratios[1, 2] = rare_var(leg_fun$row, hap_cases_pcase, maf = maf)
  ratios[1, 3] = rare_var(leg_syn$row, hap_cases_pcase, maf = maf)
  ratios[1, 4] = ratios[1, 2]/ratios[1, 3]
  ratios_case_pow = rbind(ratios_case_pow, ratios[1, 2:4])
  
  # prune the functional variants back to 100%
  rem_fun = select_var(leg_fun, exp_fun)
  hap_all_pruned = prune_var(rem_fun, hap, Nsim)
  
  # subset the pruned haplotypes (100% fun and 100% syn)
  hap_cases = sub_hap_scen(hap_all_pruned, cases, scen)
  hap_int = sub_hap_scen(hap_all_pruned, ics, scen)
  hap_cc = sub_hap_cc(hap_all_pruned, ccs)
  # subset the reference haplotypes (100% fun and 100% syn)
  hap_refs_afr = sub_hap_refs(hap_all_pruned, refs, Pop1)
  hap_refs_nfe = sub_hap_refs(hap_all_pruned, refs, Pop2)
  
  # write the haplotype files for the pruned datasets
  fwrite(hap_cases, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.cases.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_int, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.internal.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_cc, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.common.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_refs_afr, paste0(dir_out, 'chr19.block37.', Pop1, '.sim', j, '.', scen, '.ref.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  fwrite(hap_refs_nfe, paste0(dir_out, 'chr19.block37.', Pop2, '.sim', j, '.', scen, '.ref.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  ## Check ratio of fun to syn variants in cases (T1E)
  ratios[2, 2] = rare_var(leg_fun$row, hap_cases, maf = maf)
  ratios[2, 3] = rare_var(leg_syn$row, hap_cases, maf = maf)
  ratios[2, 4] = ratios[2, 2]/ratios[2, 3]
  ratios_case_t1e = rbind(ratios_case_t1e, ratios[2, 2:4])
  
  ### Check ratio of fun to syn variants in internal controls
  ratios[3, 2] = rare_var(leg_fun$row, hap_int, maf = maf)
  ratios[3, 3] = rare_var(leg_syn$row, hap_int, maf = maf)
  ratios[3, 4] = ratios[3, 2]/ratios[3, 3]
  ratios_int = rbind(ratios_int, ratios[3, 2:4])
  
  ### Check ratio of fun to syn variants in external controls (no confounding)
  ratios[4, 2] = rare_var(leg_fun$row, hap_cc, maf = maf)
  ratios[4, 3] = rare_var(leg_syn$row, hap_cc, maf = maf)
  ratios[4, 4] = ratios[4, 2]/ratios[4, 3]
  ratios_cc = rbind(ratios_cc, ratios[4, 2:4])
  
  
  #### Prune back to 80% of the functional and synonymous variants
  
  # update the allele counts for just the common controls
  leg_cc = leg
  leg_cc$count = rowSums(hap_all_pruned)
  leg_cc$MAC = ifelse(leg_cc$count>Nsim, 2*Nsim-leg_cc$count, leg_cc$count)
  
  # subset the variants
  leg_fun_cc = leg_cc %>% filter(fun=="fun")
  leg_syn_cc = leg_cc %>% filter(fun=="syn")
  
  # prune the functional and synonymous variants of the common controls to 80%
  rem_cc_fun = select_var(leg_fun_cc, exp_fun_conf)
  rem_cc_syn = select_var(leg_syn_cc, exp_syn_conf)
  hap_cc_conf = prune_var(rbind(rem_cc_fun, rem_cc_syn), hap_all_pruned, Nsim)
  
  # subset the p_conf % pruned haplotype files
  hap_cc_pruned = sub_hap_cc(hap_cc_conf, ccs)
  hap_case_pconf = sub_hap_scen(hap_cc_conf, cases, scen)
  hap_int_pruned = sub_hap_scen(hap_cc_conf, ics, scen)
  
  # write the haplotype files for the p_conf % pruned datasets
  fwrite(hap_cc_pruned, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_case_pconf, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_int_pruned, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  ### Check ratio of fun to syn variants in common controls
  ratios[5, 2] = rare_var(leg_fun$row, hap_cc_pruned, maf = maf) 
  ratios[5, 3] = rare_var(leg_syn$row, hap_cc_pruned, maf = maf) 
  ratios[5, 4] = ratios[5, 2]/ratios[5, 3]
  ratios_cc_pconf = rbind(ratios_cc_pconf, ratios[5, 2:4])
  
  ### Check ratio of fun to syn variants in cases-80%
  ratios[6, 2] = rare_var(leg_fun$row, hap_case_pconf, maf = maf) 
  ratios[6, 3] = rare_var(leg_syn$row, hap_case_pconf, maf = maf) 
  ratios[6, 4] = ratios[6, 2]/ratios[6, 3]
  ratios_case_pconf = rbind(ratios_case_pconf, ratios[6, 2:4])
  
  ### Check ratio of fun to syn variants in internal controls
  ratios[7, 2] = rare_var(leg_fun$row, hap_int_pruned, maf = maf) 
  ratios[7, 3] = rare_var(leg_syn$row, hap_int_pruned, maf = maf) 
  ratios[7, 4] = ratios[7, 2]/ratios[7, 3]
  ratios_int_pconf = rbind(ratios_int_pconf, ratios[7, 2:4])
  
  fwrite(ratios, paste0(dir_out, 'ratios_100_v_',  p_conf, '_', Ncc, '.csv'),
         quote=F, row.names=F, col.names=T, sep=',')
  
  print(j)
}
