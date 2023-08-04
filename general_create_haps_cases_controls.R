############################################################################## 
# This file is used to generate the haplotype files necessary for running
# the type I error and power calculations for proxECAT, LogProx, and iECAT-O
##############################################################################

library(data.table)
library(dplyr)

source("/home/math/siglersa/mastersProject/Input/create_haps_funcs.R")
source("/home/math/siglersa/mastersProject/Input/subset_haps_funcs.R")

Pop1 = 'AFR'
Pop2 = 'NFE'
p_case = 120
p_conf = 80
Nsim = 22500 
maf = 0.001
scen = 's1' #scenario: 's1' or 's2'

# The following have all been multiplied by 2 for number of haplotypes
afr_case_size = afr_ic_size = 8500 #s1 = 8500, s2 = 10000 
nfe_case_size = nfe_ic_size = 1500 #s1 = 1500, s2 = 0 

afr_cc_size = 17000 #10kcc = 17000, 5kcc = 8500
nfe_cc_size = 3000 #10kcc = 3000, 5kcc = 1500

afr_ref_size = nfe_ref_size = 1000

# Haplotype column indices for each pop
afr_cols = 1:38000 
nfe_cols = 38001:45000


mac_dir = '/home/math/siglersa/mastersProject/Input/'
dir_in = '/storage/math/projects/compinfo/simulations/output/NFE_AFR_pops/'
dir_out = '/home/math/siglersa/mastersProject/new_AFR_NFE_pops/cc10k/'


### read in the expected number of functional and synonymous variants from RAREsim
exp_fun = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_fun_100.txt'), header=T, sep='\t')
exp_fun_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_fun_', p_conf, '.txt'), header=T, sep='\t')
exp_syn_conf = read.table(paste0(mac_dir, 'MAC_bin_estimates_', Nsim, "_", Pop1, '_syn_', p_conf, '.txt'), header=T, sep='\t')


set.seed(1) # Will be different for each replicate but same for each run

for(j in 1:100){
  
  # read in the legend file
  leg = read.table(paste0(dir_in, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.legend'), header=T, sep='\t')
  leg$row = 1:nrow(leg)
  
  # read in the haplotype file
  hap = fread(paste0(dir_in, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.all.', p_case, 'fun.100syn.haps.gz'))
  hap = as.data.frame(hap)
  
  # add allele counts to the haplotypes
  leg$count = rowSums(hap)
  
  # subset the legend file to the functional variants (those are the only ones we'll prune)
  leg_fun = leg %>% filter(fun=="fun")
  leg_syn = leg %>% filter(fun=="syn")
  
  # Select AFR and NFE columns for necessary data sets
  cases = case_cols(afr_cols, nfe_cols, afr_case_size, nfe_case_size, scen)
  ics = int_cols(afr_cols, nfe_cols, afr_ic_size, nfe_ic_size, cases, scen)
  ccs = cc_cols(afr_cols, nfe_cols, afr_cc_size, nfe_cc_size, cases, ics, scen)
  refs = ref_cols(afr_cols, nfe_cols, afr_ref_size, nfe_ref_size, cases, ics, ccs, scen)
  
  # subset the case haplotypes (120% fun and 100% syn)
  hap_cases = sub_hap_scen(hap, cases, scen)
  
  # write the haplotype file for the cases (power)
  fwrite(hap_cases, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.cases.', p_case, 'fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # prune the functional variants back to 100%
  rem_fun = select_var(leg_fun, exp_fun)
  hap_all_pruned = prune_var(rem_fun, hap, Nsim)
  
  # subset the pruned case haplotypes (100% fun and 100% syn)
  hap_cases_pruned = sub_hap_scen(hap_all_pruned, cases, scen)
  
  # write the haplotype file for the cases (type I error)
  fwrite(hap_cases_pruned, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.cases.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the internal control haplotypes (100% fun and 100% syn)
  hap_int = sub_hap_scen(hap_all_pruned, ics, scen)
  
  # write the haplotype file for the internal controls
  fwrite(hap_int, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.internal.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the common control haplotypes (no confounding: 100% fun, 100% syn)
  hap_cc = sub_hap_cc(hap_all_pruned, ccs)
  
  #write the haplotype file for the common controls (non-confounded)
  fwrite(hap_cc, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.common.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the reference haplotypes (100% fun and 100% syn)
  hap_refs_afr = sub_hap_refs(hap_all_pruned, refs, Pop1)
  hap_refs_nfe = sub_hap_refs(hap_all_pruned, refs, Pop2)
  
  # write haplotype file for references
  fwrite(hap_refs_afr, paste0(dir_out, 'chr19.block37.', Pop1, '.sim', j, '.', scen, '.ref.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  fwrite(hap_refs_nfe, paste0(dir_out, 'chr19.block37.', Pop2, '.sim', j, '.', scen, '.ref.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  
  #### Prune back to 80% of the functional and synonymous variants
  
  # update the allele counts for just the common controls
  leg_cc = leg
  leg_cc$count = rowSums(hap_all_pruned)
  
  # subset the variants
  leg_fun_cc = leg_cc %>% filter(fun=="fun")
  leg_syn_cc = leg_cc %>% filter(fun=="syn")
  
  # prune the functional and synonymous variants of the common controls to 80%
  rem_cc_fun = select_var(leg_fun_cc, exp_fun_conf)
  rem_cc_syn = select_var(leg_syn_cc, exp_syn_conf)
  hap_cc_conf = prune_var(rbind(rem_cc_fun, rem_cc_syn), hap_all_pruned, Nsim)
  
  # subset the common controls-80%
  hap_cc_pruned = sub_hap_cc(hap_cc_conf, ccs)
  
  # write the haplotype file for the common controls-80%
  fwrite(hap_cc_pruned, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the cases-80%
  hap_case_pruned80 = sub_hap_scen(hap_cc_conf, cases, scen)
  
  # write the haplotype file for the cases-80%
  fwrite(hap_case_pruned80, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the internal controls-80%
  hap_int_pruned = sub_hap_scen(hap_cc_conf, ics, scen)
  
  # write the haplotype file for the internal controls-80%
  fwrite(hap_int_pruned, paste0(dir_out, 'chr19.block37.', Pop1, '-', Pop2, '.sim', j, '.', scen, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  print(j)
}
