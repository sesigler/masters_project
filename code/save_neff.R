# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(devtools) #need for installing Summix
library(proxecat)
library(SKAT)
library(iECAT)
# library(DescTools) # For Lin's CCC-problem installing package
# library(Summix) #need to run install_github("hendriau/Summix2")

source("/home/math/siglersa/code/functions/read_in_funcs.R")
source("/home/math/siglersa/code/functions/general_data_manip.R")
source("/home/math/siglersa/code/functions/methods_funcs.R")
source("/home/math/siglersa/code/functions/summix2_adjAF.R")
source("/home/math/siglersa/code/functions/summix2_summix.R")

# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/summix2_adjAF.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/summix2_summix.R")

# pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
# data = 'by_gene'
Pop1 = 'AFR'
Pop2 = 'NFE'
admx_pop1 = 80
admx_pop2 = 20
Nsim = '42k'
scen = 's1'
folder = '160v100v80'
p_case = 160
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 80
Ncase = Nic = 5000
Ncc = 10000  
Nref_pop1 = 704
Nref_pop2 = 642
maf = 0.001 
sim_params = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_', Pop1, 'ref', Nref_pop1, '_', Pop2, 'ref', Nref_pop2)

dir_leg = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/pruned_haps/')
dir_in = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/', sim_params, '/', scen, '/')
dir_out = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Results/Sim_', Nsim, '/', sim_params, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir_out = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Results/Sim_', Nsim, '/', sim_params, '/prox_gene_adj_', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir_out = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Results/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/', sim_params, '/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/', sim_params, '/')
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/'

# Vectors to store unadjusted p-values
neff_vec = c()


# loop through the simulation replicates
set.seed(1) 
# i=1
for (i in 1:5){
  
  # read in the legend file
  # leg = read_leg_homo(dir_leg, Pop, i)
  leg = read.table(paste0(dir_leg, 'chr19.block37.', Pop1, '_', Pop2, '.sim', i, '.', p_case, 'fun.100syn.legend'), header=T, sep='\t') #RAREsim v2.1.1 pruning only
  leg$row = 1:nrow(leg)
  
  # Need to mutate so counts get added up correctly for ZNF333
  leg = leg %>% mutate(gene = ifelse(gene == "ZNF333;ZNF333(NM_001352243:exon9:UTR5)", "ZNF333", gene))
  
  # subset the synonymous variants from the legend file
  leg_syn = leg %>% filter(fun=="syn")
  leg_fun = leg %>% filter(fun=="fun")
  
  # read in the haplotype files
  hap_case = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn)
  hap_ic = read_hap(dir_in, Pop1, Pop2, i, scen, "internal.controls", p_int_fun, p_int_syn)
  hap_cc = read_hap(dir_in, Pop1, Pop2, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  hap_ref_pop1 = read_ref(dir_in, Pop1, i, scen, p_fun = 100, p_syn = 100)
  hap_ref_pop2 = read_ref(dir_in, Pop2, i, scen, p_fun = 100, p_syn = 100)
  
  # convert the haplotypes into genotypes
  geno_case = make_geno(hap_case)
  geno_ic = make_geno(hap_ic)
  geno_cc = make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_case = calc_allele_freqs(geno_case, Ncase, Pop=NULL)
  count_ic = calc_allele_freqs(geno_ic, Nic, Pop=NULL)
  count_cc = calc_allele_freqs(geno_cc, Ncc, Pop=NULL)
  
  count_ref_pop1 = calc_allele_freqs(hap_ref_pop1, Nref_pop1, Pop=Pop1)
  count_ref_pop2 = calc_allele_freqs(hap_ref_pop2, Nref_pop2, Pop=Pop2)
  
  # Commbine data with references for Summix
  cc_refs = cbind(count_cc, count_ref_pop1, count_ref_pop2)
  case_refs = cbind(count_case, count_ref_pop1, count_ref_pop2)
  
  # Estimate ancestry proportions using only COMMON variants
  cc_est_prop = est_props(cc_refs, Pop1, Pop2, maf)
  case_est_prop = est_props(case_refs, Pop1, Pop2, maf)
  
  # Add row index column to cc_refs in case summix removes variants during adjustment 
  cc_refs$row <- 1:nrow(cc_refs)
  
  # Calculate adjusted AFs
  adj_Neff = calc_adjusted_AF(cc_refs, Pop1, Pop2, case_est_prop, cc_est_prop, Nref=c(Nref_pop1, Nref_pop2), Ncc, Neff=TRUE)
  
  # return effective sample size from Neff adjusted data
  Neff = adj_Neff[[2]]
  
  neff_vec = c(neff_vec, Neff)
  
  
  print(i)
}


# Save the effective sample sizes
write.csv(neff_vec, paste0(dir_out, "neff.csv"), quote=F, row.names=F)

