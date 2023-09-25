############################################################################## 
# This file is used to generate the expected MAC bin estimates for a specified
# population and different levels of confounding, but not stratified by 
# functional or synonymous status
##############################################################################

library(RAREsim)

pop <- 'NFE'
nsim <- 20000
reg_size <- 19.029
p_case <- 120
p_exp <- 100
p_conf <- 80

dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/'


### total number of variants in the region
exp_var <- reg_size*nvariant(N = nsim, pop = pop)

### minor allele frequencies
maf1 <- round(0.01*(2*nsim))
maf0.5 <- round(maf1/2)
maf0.25 <- round(maf0.5/2)

### minor allele count bins
if (nsim > 3500){
  mac <- data.frame(Lower = c(1, 2, 3, 6, 11, 21, maf0.5+1),
                    Upper = c(1, 2, 5, 10, 20, maf0.5, maf1))
} else {
  mac <- data.frame(Lower = c(1, 2, 3, 6, maf0.25+1, maf0.5+1),
                    Upper = c(1, 2, 5, maf0.25, maf0.5, maf1))
}

### allele frequency spectrum function
bin_props <- afs(mac_bins = mac, pop = pop)

### expected number of variants per bin
bin_estimates_case <- expected_variants(Total_num_var = (p_case/100)*exp_var, mac_bin_prop = bin_props)
bin_estimates_exp <- expected_variants(Total_num_var = exp_var, mac_bin_prop = bin_props)
bin_estimates_conf <- expected_variants(Total_num_var = (p_conf/100)*exp_var, mac_bin_prop = bin_props)

### write file with variants to prune by removing
write.table(bin_estimates_case, paste0(dir_out, 'MAC_bin_estimates_', nsim, '_', pop, '_', p_case, '.txt'), row.names=F, quote=F, sep = '\t')
write.table(bin_estimates_exp, paste0(dir_out, 'MAC_bin_estimates_', nsim, '_', pop, '_', p_exp, '.txt'), row.names=F, quote=F, sep = '\t')
write.table(bin_estimates_conf, paste0(dir_out, 'MAC_bin_estimates_', nsim, '_', pop, '_', p_conf, '.txt'), row.names=F, quote=F, sep = '\t')
