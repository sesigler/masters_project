############################################################################## 
# This file is used to generate the expected MAC bin estimates for a specified
# population and different levels of confounding
##############################################################################

library(RAREsim)
library(dplyr)

Pop = 'NFE'
p_case = 120
p_conf = 99
Nsim = 20000
reg_size = 19.029 


dir_leg = paste0("/storage/math/projects/RAREsim/Cases/Sim_20k/", Pop, "/data/")
dir = "/storage/math/projects/compinfo/simulations/input/"

dir_leg = "C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/"
dir = "C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/input/"

### prep the legend file
leg = read.table(paste0(dir_leg, "chr19.block37.", Pop, ".sim1.copy.legend"), sep='\t', header=T)
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

### write file with bin estimates
write.table(bins_fun_cases, paste0(dir, "MAC_bin_estimates_", Nsim, "_", Pop, "_fun_", p_case, ".txt"), row.names=F, quote=F, sep='\t')
write.table(bins_fun_exp, paste0(dir, "MAC_bin_estimates_", Nsim, "_", Pop, "_fun_100.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_fun_conf, paste0(dir, "MAC_bin_estimates_", Nsim, "_", Pop, "_fun_", p_conf, ".txt"), row.names=F, quote=F, sep='\t')
write.table(bins_syn_exp, paste0(dir, "MAC_bin_estimates_", Nsim, "_", Pop, "_syn_100.txt"), row.names=F, quote=F, sep='\t')
write.table(bins_syn_conf, paste0(dir, "MAC_bin_estimates_", Nsim, "_", Pop, "_syn_", p_conf, ".txt"), row.names=F, quote=F, sep='\t')
