############################################################################## 
# This file is used to generate the expected MAC bin estimates for a specified
# population and different levels of confounding
# This file also combines the last two mac bins to avoid the error when the 
# observed is less than the expected in the last mac bin
##############################################################################

#' combine_bins
#' 
#' @description
#' Combine the last two rows of the mac bin estimates
#' 
#' @param bins the estimated mac bins
#' 
#' @return the new mac bins with the 6th and 7th bins added together

combine_bins = function(bins) {
  
  # Make a copy of the mac bins but only keep the first 5 bins
  bins2 = bins[1:5,]
  
  # Keep Lower same as 6th bin, Upper value same as 7th bin, and add Expected_var values from bins 6 and 7
  bins2[6, 1:3] <- c(bins[6, 1], bins[7, 2], bins[6, 3] + bins[7, 3]) 
  
  return(bins2)
}

library(RAREsim)
library(dplyr)

Pop = 'AFR'
Nsim_hapgen = "100000"
haps_fold = 100
p_case = 160
p_conf = 80
Nsim = 42000
reg_size = 19.029 


# dir_leg = paste0("/storage/math/projects/RAREsim/Cases/Sim_20k/", Pop, "/data/")
dir_leg = paste0('/home/math/siglersa/Sim_', haps_fold, 'k/', Pop, '/data/')
dir = '/storage/math/projects/compinfo/simulations/input/'
dir_out = paste0('/home/math/siglersa/mac_bin_estimates/', Pop, '/')

# dir_leg = "C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/mac_bin_data/"
# dir = "C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/mac_bin_data/"

### prep the legend file
leg = read.table(paste0(dir_leg, 'chr19.block37.', Pop, '.sim1.', Nsim_hapgen, '.copy.legend'), sep='\t', header=T)
#leg$fun = ifelse(leg$fun=="synonymous SNV", "syn", "fun")
#leg$exonic[grepl("exonic", leg$exonic)] = "exonic"
#leg$gene[grepl("ZNF333", leg$gene)] = "ZNF333"

### number of variants target data
nvar = read.table(paste0(dir, "Block37_fun_syn_num_var.txt"), header=T) %>% 
  filter(pop==tolower(Pop)) %>% rename(n=downsample)

### divide by the region size for per_kb
nvar2 = nvar %>% mutate(fun_per_kb=obs_fun/reg_size, syn_per_kb=obs_syn/reg_size)

### estimate parameters by functional status for nvariant function
nvar_fun = fit_nvariant(nvar2 %>% select(n, fun_per_kb))
nvar_syn = fit_nvariant(nvar2 %>% select(n, syn_per_kb))

### total number of variants in the region
exp_var_fun = reg_size*nvariant(phi = nvar_fun$phi, omega = nvar_fun$omega, N = Nsim)
exp_var_syn = reg_size*nvariant(phi = nvar_syn$phi, omega = nvar_syn$omega, N = Nsim)

### minor allele frequencies for the target data
Ntar = last(nvar$n)
tar_maf1 = round(0.01*(2*Ntar))
tar_maf0.5 = round(tar_maf1/2)
tar_maf0.25 = round(tar_maf0.5/2)

### minor allele count bins for the target data
if (Ntar > 3500){
  mac_tar = data.frame(Lower = c(1, 2, 3, 6, 11, 21, tar_maf0.5+1),
                       Upper = c(1, 2, 5, 10, 20, tar_maf0.5, tar_maf1))
} else {
  mac_tar = data.frame(Lower = c(1, 2, 3, 6, tar_maf0.25+1, tar_maf0.5+1),
                       Upper = c(1, 2, 5, tar_maf0.25, tar_maf0.5, tar_maf1))
}

### subset legend file to variants observed in gnomad (target data)
leg2 = leg %>% filter(AC!='.', AC!=0, exonic=="exonic") # Megan just filtered by exonic
leg2$AC = as.numeric(leg2$AC)

### subset target data by functional status
leg_fun = leg2 %>% filter(fun=="fun")
leg_syn = leg2 %>% filter(fun=="syn")

### count the number of variants within each MAC bin (target data)
mac_fun_tar = mac_syn_tar = mac_tar
mac_fun_tar$count = mac_syn_tar$count = 0
for (k in 1:nrow(mac_tar)){
  mac_fun_tar$count[k] = sum(between(leg_fun$AC, mac_tar$Lower[k], mac_tar$Upper[k]))
  mac_syn_tar$count[k] = sum(between(leg_syn$AC, mac_tar$Lower[k], mac_tar$Upper[k]))
}

### convert the counts to proportions
mac_fun_tar2 = mac_fun_tar %>% mutate(Prop = count/nrow(leg_fun)) %>% select(-count)
mac_syn_tar2 = mac_syn_tar %>% mutate(Prop = count/nrow(leg_syn)) %>% select(-count)

#write.table(mac_fun_tar2, 'AFS_fun_target_data.txt', row.names = FALSE, quote = FALSE, sep = '\t')
#write.table(mac_syn_tar2, 'AFS_syn_target_data.txt', row.names = FALSE, quote = FALSE, sep = '\t')

### allele frequency spectrum function
afs_fun = fit_afs(mac_fun_tar2)
afs_syn = fit_afs(mac_syn_tar2)

### minor allele frequencies for the simulated data
sim_maf1 = round(0.01*(2*Nsim))
sim_maf0.5 = round(sim_maf1/2)
sim_maf0.25 = round(sim_maf0.5/2)

### minor allele count bins for the simulated data
if (Nsim > 3500){
  mac_sim = data.frame(Lower = c(1, 2, 3, 6, 11, 21, sim_maf0.5+1),
                       Upper = c(1, 2, 5, 10, 20, sim_maf0.5, sim_maf1))
} else {
  mac_sim = data.frame(Lower = c(1, 2, 3, 6, sim_maf0.25+1, sim_maf0.5+1),
                       Upper = c(1, 2, 5, sim_maf0.25, sim_maf0.5, sim_maf1))
}
mac_fun_sim = mac_syn_sim = mac_sim

### add the proportions from the target data
mac_fun_sim$Prop = afs_fun$Fitted_results$Prop
mac_syn_sim$Prop = afs_syn$Fitted_results$Prop

### expected number of variants per bin
bins_fun_cases = expected_variants(Total_num_var = (p_case/100)*exp_var_fun, mac_bin_prop = mac_fun_sim)
bins_fun_exp = expected_variants(Total_num_var = exp_var_fun, mac_bin_prop = mac_fun_sim)
bins_fun_conf = expected_variants(Total_num_var = (p_conf/100)*exp_var_fun, mac_bin_prop = mac_fun_sim)
bins_syn_exp = expected_variants(Total_num_var = exp_var_syn, mac_bin_prop = mac_syn_sim)
bins_syn_conf = expected_variants(Total_num_var = (p_conf/100)*exp_var_syn, mac_bin_prop = mac_syn_sim)

# Calculate combined mac bin estimates
bins_fun_cases2 <- combine_bins(bins_fun_cases)
bins_fun_exp2 <- combine_bins(bins_fun_exp)
bins_fun_conf2 <- combine_bins(bins_fun_conf)
bins_syn_exp2 <- combine_bins(bins_syn_exp)
bins_syn_conf2 <- combine_bins(bins_syn_conf)

### write file with bin estimates
write.table(bins_fun_cases, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_fun_", p_case, ".txt"), row.names=F, quote=F, sep='\t')
write.table(bins_fun_exp, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_fun_100.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_fun_conf, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_fun_", p_conf, ".txt"), row.names=F, quote=F, sep='\t')
write.table(bins_syn_exp, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_syn_100.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_syn_conf, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_syn_", p_conf, ".txt"), row.names=F, quote=F, sep='\t')

### write file with 6 bin estimates
write.table(bins_fun_cases2, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_fun_", p_case, "_6bins.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_fun_exp2, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_fun_100_6bins.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_fun_conf2, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_fun_", p_conf, "_6bins.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_syn_exp2, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_syn_100_6bins.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_syn_conf2, paste0(dir_out, "MAC_bin_estimates_", format(Nsim, scientific=F), "_", Pop, "_syn_", p_conf, "_6bins.txt"), row.names=F, quote=F, sep='\t')

