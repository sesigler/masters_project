# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(proxecat)
# library(devtools) #need for installing Summix
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
scen = 's2'
folder = '160v100v80'
p_case = 160
p_case_fun = p_case_syn = p_int_fun = p_int_syn = p_exp = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 80
Ncase = Nic = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
Nref = 10000
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
# pi_tar1 = 1 #pi.target for AFR: 0.80 for s1 or 1 for s2
# pi_tar2 = 0 #pi.target for NFE: 0.20 for s1 or 0 for s2
# genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)


dir_leg = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Sim_42k/', folder, '/')
dir_in = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Sim_42k/', folder, '/datasets/', scen, '/')
dir_out = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Results/Sim_42k/prox_gene_adj_', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir_out = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Results/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/')
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/'

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/admixed/'
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')

# Vectors to store p-values
# proxECAT
prox_int_p = prox_ext_p = prox_ext_p_var_adj = prox_ext_p_gene_adj = c() 

# proxECAT-weighted
proxW_int_p = proxW_ext_p = proxW_ext_p_var_adj = proxW_ext_p_gene_adj = c() 



# loop through the simulation replicates
set.seed(1) 
# i=1
for (i in 1:5){
  
  # read in the legend file
  # leg = read_leg_homo(dir_leg, Pop, i)
  leg = read.table(paste0(dir_leg, 'chr19.block37.', Pop1, '_', Pop2, '.sim', i, '.', p_case, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t') #RAREsim v2.1.1 pruning only
  leg$row = 1:nrow(leg)
  
  # Need to mutate so counts get added up correctly for ZNF333
  leg = leg %>% mutate(gene = ifelse(gene == "ZNF333;ZNF333(NM_001352243:exon9:UTR5)", "ZNF333", gene))
  
  # subset the synonymous variants from the legend file
  leg_syn = leg %>% filter(fun=="syn")
  leg_fun = leg %>% filter(fun=="fun")
  
  # read in the haplotype files
  hap_case = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn)
  # hap_cases_power = read_hap_homo(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn) # pcase % fun 100% syn
  # hap_cases_t1e = read_hap_homo(dir_in, Pop1, Pop2, i, scen, "cases", p_int_fun, p_int_syn) # 100% fun 100% syn
  hap_ic = read_hap(dir_in, Pop1, Pop2, i, scen, "internal.controls", p_int_fun, p_int_syn)
  hap_cc = read_hap(dir_in, Pop1, Pop2, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  hap_ref_pop1 = read_ref(dir_in, Pop1, i, scen, p_exp, p_exp)
  hap_ref_pop2 = read_ref(dir_in, Pop2, i, scen, p_exp, p_exp)
  
  # FOR POWER ONLY
  # Create a new hap cases dataframe that merges the cases used for power and t1e but only contains
  # the genes associated with each calculation
  # hap_cases = merge_cases(hap_cases_power, hap_cases_t1e, leg, genes_power)
  
  # convert the haplotypes into genotypes
  geno_case = make_geno(hap_case)
  geno_ic = make_geno(hap_ic)
  geno_cc = make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_case = calc_allele_freqs(geno_case, Ncase, Pop=NULL)
  count_ic = calc_allele_freqs(geno_ic, Nic, Pop=NULL)
  count_cc = calc_allele_freqs(geno_cc, Ncc, Pop=NULL)
  
  count_ref_pop1 = calc_allele_freqs(hap_ref_pop1, Nref, Pop=Pop1)
  count_ref_pop2 = calc_allele_freqs(hap_ref_pop2, Nref, Pop=Pop2)
  
  # Commbine data with references for Summix
  cc_refs = cbind(count_cc, count_ref_pop1, count_ref_pop2)
  case_refs = cbind(count_case, count_ref_pop1, count_ref_pop2)
  
  # Estimate ancestry proportions using only COMMON variants
  cc_est_prop = est_props(cc_refs, Pop1, Pop2, maf)
  case_est_prop = est_props(case_refs, Pop1, Pop2, maf)
  # int_est_prop = est_props(int_refs, Pop1, Pop2, maf)
  
  # prop_ests_cc <- rbind(prop_ests_cc, cc_est_prop)
  # prop_ests_cases <- rbind(prop_ests_cases, cases_est_prop)
  # prop_ests_int <- rbind(prop_ests_int, int_est_prop)
  
  # Add row index column to cc_refs in case summix removes variants during adjustment 
  cc_refs$row <- 1:nrow(cc_refs)
  
  # Calculate adjusted AFs
  count_cc_adj = calc_adjusted_AF(cc_refs, Pop1, Pop2, case_est_prop, cc_est_prop, Nref, Ncc)
  
  # Identify variants where AF >= 1-maf
  flip_int = which(count_case$af >= 1-maf | count_ic$af >= 1-maf)
  flip_ext = which(count_case$af >= 1-maf | count_cc$af >= 1-maf)
  flip_ext_adj = which(count_case$af >= 1-maf | count_cc_adj$af >= 1-maf)
  
  # Flip the data at the variants identified above for all the different combination of datasets
  # If no variants need to be flipped, return the unaltered datasets
  # Cases and internal controls
  int_data = flip_data(leg, flip_int, geno_case, count_case, Ncase, cntrl="int", geno_ic, count_ic, Nic, 
                       geno.cc=NULL, count.cc=NULL, count.cc.adj=NULL, Ncc=NULL, adj=FALSE)
  
  leg_int = int_data[[1]]
  geno_case_int = int_data[[2]]
  geno_ic_int = int_data[[3]]
  count_case_int = int_data[[5]]
  count_ic_int = int_data[[6]]
  
  # Cases and external controls
  ext_data = flip_data(leg, flip_ext, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                       geno_cc, count_cc, count.cc.adj=NULL, Ncc, adj=FALSE)
  
  leg_ext = ext_data[[1]]
  geno_case_ext = ext_data[[2]]
  geno_cc_ext = ext_data[[4]]
  count_case_ext = ext_data[[5]]
  count_cc_ext = ext_data[[7]]
  
  # Cases and adjusted external controls
  ext_adj_data = flip_data(leg, flip_ext_adj, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                           geno.cc=NULL, count.cc=NULL, count_cc_adj, Ncc, adj=TRUE)
  
  leg_ext_adj = ext_adj_data[[1]]
  geno_case_ext_adj = ext_adj_data[[2]]
  count_case_ext_adj = ext_adj_data[[5]]
  count_cc_ext_adj = ext_adj_data[[8]]
  
  # identify the common variants
  common_int = leg[which(count_case_int$af > maf | count_ic_int$af > maf),]
  common_ext = leg[which(count_case_ext$af > maf | count_cc_ext$af > maf),]
  common_ext_adj = leg[which(count_case_ext_adj$af > maf | count_cc_ext_adj$af > maf),]
  
  # proxECAT
  counts_int_wide = prox_gene_data_prep(count_case_int, count_ic_int, leg_int, common_int)
  counts_ext_wide = prox_gene_data_prep(count_case_ext, count_cc_ext, leg_ext, common_ext)
  counts_ext_wide_adj = prox_gene_data_prep(count_case_ext_adj, count_cc_ext_adj, leg_ext_adj, common_ext_adj)
  
  # Store the proxECAT and proxECAT-weighted p-values
  prox_int_p = rbind(prox_int_p, counts_int_wide$prox)
  prox_ext_p = rbind(prox_ext_p, counts_ext_wide$prox)
  prox_ext_p_var_adj = rbind(prox_ext_p_var_adj, counts_ext_wide_adj$prox)
  
  proxW_int_p = rbind(proxW_int_p, counts_int_wide$prox_w)
  proxW_ext_p = rbind(proxW_ext_p, counts_ext_wide$prox_w)
  proxW_ext_p_var_adj = rbind(proxW_ext_p_var_adj, counts_ext_wide_adj$prox_w)
  
  # Adjust AFs at the gene level instead of variant level
  # Combine case and control data
  count_ext_ref <- cbind(count_case, cc_refs)
  colnames(count_ext_ref)[1:4] <- c("ac_case", "af_case", "ac_cc", "af_cc")
  count_ext_ref2 = count_ext_ref %>% mutate(gene = leg$gene, id = leg$id, fun = leg$fun)

  # Identify common variants
  common_ext2 = leg[which((count_ext_ref2$af_case > maf & count_ext_ref2$af_case < 1-maf) | (count_ext_ref2$af_cc > maf & count_ext_ref2$af_cc < 1-maf)),]

  # Filter out common variants
  count_ext_rare = count_ext_ref2 %>% filter(!(id %in% common_ext2$id)) %>% mutate(across(all_of(c("gene","id", "fun")), as.factor))

  # Subset data by functional status
  data_ext_syn = count_ext_rare %>% filter(fun == "syn")
  data_ext_fun = count_ext_rare %>% filter(fun == "fun")

  # Only needed if dividng by # of fun/syn variants per gene
  # syn_gene = data.frame(table(data_ext_syn$gene))
  # fun_gene = data.frame(table(data_ext_fun$gene))

  # Sum up the ACs and AFs
  data_ext_syn2 = data_ext_syn %>% group_by(gene) %>% summarise(af_case = sum(af_case), af_cc = sum(af_cc), af_afr = sum(af_afr), af_nfe = sum(af_nfe),
                                                                ac_case = sum(ac_case), ac_cc = sum(ac_cc), ac_afr = sum(ac_afr), ac_nfe = sum(ac_nfe)) 

  # Perform the adjustment by gene
  adj_AF_syn <- adjAF(data = data_ext_syn2,
                      reference = c("af_afr", "af_nfe"),
                      observed = "af_cc",
                      pi.target = c(case_est_prop[, "af_afr"], case_est_prop[, "af_nfe"]),
                      pi.observed = c(cc_est_prop[, "af_afr"], cc_est_prop[, "af_nfe"]),
                      adj_method = "average",
                      N_reference = c(Nref, Nref),
                      N_observed = Ncc,
                      filter = TRUE)

  # Add adj AF to data frame
  data_ext_syn2$adj_af <- adj_AF_syn$adjusted.AF$adjustedAF

  # Calculate the adjusted AC
  data_ext_syn2$adj_ac <- round(data_ext_syn2$adj_af*(2*Ncc))

  # Sum up ACs and AFs for fun variants
  data_ext_fun2 = data_ext_fun %>% group_by(gene) %>% summarise(af_case = sum(af_case), af_cc = sum(af_cc), af_afr = sum(af_afr), af_nfe = sum(af_nfe),
                                                                ac_case = sum(ac_case), ac_cc = sum(ac_cc), ac_afr = sum(ac_afr), ac_nfe = sum(ac_nfe))

  # Adjust the fun common controls by gene
  adj_AF_fun <- adjAF(data = data_ext_fun2,
                      reference = c("af_afr", "af_nfe"),
                      observed = "af_cc",
                      pi.target = c(case_est_prop[, "af_afr"], case_est_prop[, "af_nfe"]),
                      pi.observed = c(cc_est_prop[, "af_afr"], cc_est_prop[, "af_nfe"]),
                      adj_method = "average",
                      N_reference = c(Nref, Nref),
                      N_observed = Ncc,
                      filter = TRUE)

  # Add adj AF to data frame
  data_ext_fun2$adj_af <- adj_AF_fun$adjusted.AF$adjustedAF

  # Calculate the adjusted AC
  data_ext_fun2$adj_ac <- round(data_ext_fun2$adj_af*(2*Ncc))

  # Combine fun and syn data
  colnames(data_ext_syn2) <- paste(colnames(data_ext_syn2), "syn", sep = "_")
  colnames(data_ext_fun2) <- paste(colnames(data_ext_fun2), "fun", sep = "_")
  data_prox_adj <- cbind(data_ext_fun2[, c("gene_fun", "ac_case_fun", "adj_ac_fun")], data_ext_syn2[, c("ac_case_syn", "adj_ac_syn")])

  # Run proxECAT and proxECAT-weighted
  counts_gene_adj = data_prox_adj %>% mutate(case_ratio = ac_case_fun/ac_case_syn,
                                             control_ratio = adj_ac_fun/adj_ac_syn,
                                             case_fun_w = ac_case_fun/median(case_ratio),
                                             control_fun_w = adj_ac_fun/median(control_ratio)) %>%
    mutate(prox = ifelse((ac_case_fun + adj_ac_fun < 5) | (ac_case_syn + adj_ac_syn < 5), NA,
                         proxecat(ac_case_fun, ac_case_syn, adj_ac_fun, adj_ac_syn)$p.value),
           prox_w = ifelse((case_fun_w + control_fun_w < 5) | (ac_case_syn + adj_ac_syn < 5), NA,
                           proxecat(case_fun_w, ac_case_syn, control_fun_w, adj_ac_syn)$p.value))
  
  # Store the proxECAT and proxECAT-weighted p-values
  prox_ext_p_gene_adj = rbind(prox_ext_p_gene_adj, counts_gene_adj$prox)
  proxW_ext_p_gene_adj = rbind(proxW_ext_p_gene_adj, counts_gene_adj$prox_w)
  
  
  print(i)
}

# Get gene names
genes = levels(droplevels(as.factor(leg$gene)))

# Set col names to the genes
colnames(prox_int_p) = colnames(prox_ext_p) = colnames(prox_ext_p_var_adj) = colnames(prox_ext_p_gene_adj) = genes
colnames(proxW_int_p) = colnames(proxW_ext_p) = colnames(proxW_ext_p_var_adj) = colnames(proxW_ext_p_gene_adj) = genes

# Set file path name
file_path = paste0(int_prune, "_v_", ext_prune, "_", Pop1, "_", Pop2, "_", scen, "_maf", maf, ".txt")

# ProxECAT
write.table(prox_int_p, paste0(dir_out, "T1e_gene_prox_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_p, paste0(dir_out, "T1e_gene_prox_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_p_var_adj, paste0(dir_out, "T1e_gene_prox_ext_var_adj_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_p_gene_adj, paste0(dir_out, "T1e_gene_prox_ext_gene_adj_", file_path), quote=F, row.names=F, col.names=T)
# ProxECAT-weighted
write.table(proxW_int_p, paste0(dir_out, "T1e_gene_prox_weighted_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_p, paste0(dir_out, "T1e_gene_prox_weighted_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_p_var_adj, paste0(dir_out, "T1e_gene_prox_weighted_ext_var_adj_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_p_gene_adj, paste0(dir_out, "T1e_gene_prox_weighted_ext_gene_adj_", file_path), quote=F, row.names=F, col.names=T)

# library(ggplot2)
# g1 <- ggplot(data_ext2 %>% filter(fun == "fun"), aes(x=adj_maf, y=af_case)) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1) +
#   theme_bw() +
#   xlab("Adjusted AF") + ylab("AFR AF") +
#   geom_text(mapping = aes(x = 0.002, y = 0.007),
#             label = paste0("CCC fun = ", round(DescTools::CCC(data_ext_fun2$adj_maf, data_ext_fun2$af_case)$rho.c$est, 6)))
# g1
