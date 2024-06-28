################################################################################ 
# This file is used to create the datasets necessary for running type I error 
# and power calculations for proxECAT, LogProx, and iECAT-O on an ADMIXED 
# population consisting of any number of continental ancestry populations
################################################################################
# Note: The order of admx_props_int, admx_props_ext, Nrefs, and Nhaps should all
# match up with the order of Pops
################################################################################

library(data.table)
library(dplyr)

source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/create_haps_funcs.R")

# Set simulation parameters
Pop_admx = 'LTX'
Pops = c('IAM', 'NFE', 'EAS', 'AFR')
admx_props_int= c(47, 44, 5, 4) # internal sample admixture proportions
admx_props_ext= c(47, 44, 5, 4) # common control admixture proportions
p_case = 160
p_conf = 80
scen = 's1'
sub_scen = 'default'
Ncase = 2000
Nic = 2000
Ncc = 10000
Nrefs = c(2000, 2000, 2000, 2000)
Nhaps = c(17160, 16320, 5400, 5120)

dir_in = paste0('/home/math/siglersa/admixed/', paste(paste(admx_props_ext, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/pruned_haps/')
dir_out = paste0('/home/math/siglersa/admixed/', paste(paste(admx_props_ext, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/datasets/')

# dir_in = dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/', scen, '/', sub_scen, '/')
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_', Nsim, '/pruned_haps/')

# Create empty list to store number of haplotypes for each dataset by pop
Ncase_pops <- setNames(vector("list", length(Pops)), paste0("Ncase_pop", 1:length(Pops)))
Nic_pops <- setNames(vector("list", length(Pops)), paste0("Nic_pop", 1:length(Pops)))
Ncc_pops <- setNames(vector("list", length(Pops)), paste0("Ncc_pop", 1:length(Pops)))
Nref_pops <- setNames(vector("list", length(Pops)), paste0("Nref_pop", 1:length(Pops)))

# Calculate and assign number of haplotypes for each dataset by population
for (i in 1:length(Pops)) {
  Ncase_pops[i] <- Ncase*2*(admx_props_int[i]/100)
  Nic_pops[i] <- Nic*2*(admx_props_int[i]/100)
  Ncc_pops[i] <- Ncc*2*(admx_props_ext[i]/100)
  Nref_pops[i] <- Nrefs[i]*2
}

# Create empty list to store column numbers corresponding to each pop
pop_cols <- setNames(vector("list", length(Pops)), paste0("pop", 1:length(Pops), "_cols"))

# Assign the column numbers for the first pop
pop_cols[[1]] <- 1:Nhaps[1]

# Assign column numbers for remaining pops if there's more than one pop
if (length(Pops) > 1) {
  # Start at 2 because we already assigned the column numbers for the first pop
  for (i in 2:length(pop_cols)) {
    
    # Starting column index of pop i
    start <- sum(Nhaps[1:(i-1)]) + 1 
    
    # ending column index for pop i
    end <- sum(Nhaps[1:i]) 
    
    # save column positions for pop i
    pop_cols[[i]] <- start:end 
  }
}


set.seed(1) # Will be different for each replicate but same for each run
# j = 1
for(j in 1:100){
  
  # Read in legend files
  leg_pcase = read.table(paste0(dir_in, 'chr19.block37.', Pop_admx, '.sim', j, '.', p_case, 'fun.100syn.legend'), header=T, sep='\t')
  leg_pcase$row = 1:nrow(leg_pcase)
  
  leg_pexp = read.table(paste0(dir_in, 'chr19.block37.', Pop_admx, '.sim', j, '.100fun.100syn.legend'), header=T, sep='\t')
  
  leg_pconf = read.table(paste0(dir_in, 'chr19.block37.', Pop_admx, '.sim', j, '.', p_conf, 'fun.', p_conf, 'syn.legend'), header=T, sep='\t')
  
  ### Read in pruned hap files
  hap_pcase = fread(paste0(dir_in, 'chr19.block37.', Pop_admx, '.sim', j, '.all.', p_case, 'fun.100syn.haps.gz'))
  hap_pcase = as.data.frame(hap_pcase)
  
  hap_pexp = fread(paste0(dir_in, 'chr19.block37.', Pop_admx, '.sim', j, '.all.100fun.100syn.haps.gz'))
  hap_pexp = as.data.frame(hap_pexp)
  
  hap_pconf = fread(paste0(dir_in, 'chr19.block37.', Pop_admx, '.sim', j, '.all.', p_conf, 'fun.', p_conf, 'syn.haps.gz'))
  hap_pconf = as.data.frame(hap_pconf)
  
  # Add rows of zeros back into 100% and p_conf % pruned hap files
  hap_exp_pruned = data.frame(matrix(0, nrow=nrow(leg_pcase), ncol=ncol(hap_pcase)))
  hap_exp_pruned[which(leg_pcase$id %in% leg_pexp$id),] = hap_pexp
  
  hap_all_pruned = data.frame(matrix(0, nrow=nrow(leg_pcase), ncol=ncol(hap_pcase)))
  hap_all_pruned[which(leg_pcase$id %in% leg_pconf$id),] = hap_pconf
  
  # Randomly select the columns (i.e. individuals) for each dataset
  
  # Create empty list to store columns for each dataset
  case_cols <- setNames(vector("list", length(Pops)), paste0("case_cols_pop", 1:length(Pops)))
  ic_cols <- setNames(vector("list", length(Pops)), paste0("ic_cols_pop", 1:length(Pops)))
  cc_cols <- setNames(vector("list", length(Pops)), paste0("cc_cols_pop", 1:length(Pops)))
  ref_cols <- setNames(vector("list", length(Pops)), paste0("ref_cols_pop", 1:length(Pops)))
  
  # Randomly sample the columns for each pop in each dataset
  for (i in 1:length(Pops)) {
    case_cols[[i]] <- sort(sample(x=pop_cols[[i]], size = Ncase_pops[[i]], replace = FALSE))
    ic_cols[[i]] <- sort(sample(x=pop_cols[[i]][! pop_cols[[i]] %in% case_cols[[i]]], size = Nic_pops[[i]], replace = FALSE))
    cc_cols[[i]] <- sort(sample(x=pop_cols[[i]][! pop_cols[[i]] %in% c(case_cols[[i]], ic_cols[[i]])], size = Ncc_pops[[i]], replace = FALSE))
    ref_cols[[i]] <- sort(sample(x=pop_cols[[i]][! pop_cols[[i]] %in% c(case_cols[[i]], ic_cols[[i]], cc_cols[[i]])], size = Nref_pops[[i]], replace = FALSE))
  }
  
  # subset the case haplotypes (pcase % fun and 100% syn)
  hap_cases_pcase = hap_pcase[, c(unlist(case_cols, use.names = FALSE))]
  
  # write the haplotype file for the cases (power)
  fwrite(hap_cases_pcase, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.cases.', p_case, 'fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # subset the pruned haplotypes for 100% pruned
  hap_cases = hap_exp_pruned[, c(unlist(case_cols, use.names = FALSE))]
  hap_int = hap_exp_pruned[, c(unlist(ic_cols, use.names = FALSE))]
  hap_cc = hap_exp_pruned[, c(unlist(cc_cols, use.names = FALSE))]
  
  # write the haplotype files for the cases, internal and common controls
  fwrite(hap_cases, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.cases.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_int, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.internal.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_cc, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.common.controls.100fun.100syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  # Create list to store reference haplotypes
  ref_haps <- setNames(vector("list", length(Pops)), paste0("hap_refs_pop", 1:length(Pops)))
  
  for (i in 1:length(Pops)) {
    
    # Subset the pruned ref haps for 100% pruned
    ref_haps[[i]] <- hap_exp_pruned[, ref_cols[[i]]] #CHECK THIS
    
    # Save the ref hap file for each pop
    fwrite(ref_haps[[i]], paste0(dir_out, 'chr19.block37.', Pops[i], '.sim', j, '.', scen, '.refs.100fun.100syn.haps.gz'),
           quote=F, row.names=F, col.names=F, sep=' ')
  }
  
  # Subset the datasets for the p_conf % pruned haplotype
  hap_cases_pconf = hap_all_pruned[, c(unlist(case_cols, use.names = FALSE))]
  hap_int_pconf = hap_all_pruned[, c(unlist(ic_cols, use.names = FALSE))]
  hap_cc_pconf = hap_all_pruned[, c(unlist(cc_cols, use.names = FALSE))]
  
  # write the haplotype files for the cases, internal and common controls p_conf % pruned
  fwrite(hap_cases_pconf, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.cases.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_int_pconf, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.internal.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  fwrite(hap_cc_pconf, paste0(dir_out, 'chr19.block37.', Pop_admx, '.sim', j, '.', scen, '.common.controls.', p_conf, 'fun.', p_conf, 'syn.haps.gz'),
         quote=F, row.names=F, col.names=F, sep=' ')
  
  print(j)
}