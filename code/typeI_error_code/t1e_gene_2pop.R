# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(devtools) #need for installing Summix
library(proxecat)
library(SKAT)
library(iECAT)
library(Summix) #need to run install_github("hendriau/Summix2")

source("/home/math/siglersa/mastersProject/Input/read_in_funcs.R")
source("/home/math/siglersa/mastersProject/Input/general_data_manip.R")
source("/home/math/siglersa/mastersProject/Input/methods_funcs.R")

source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")

Pop1 = 'AFR'
Pop2 = 'NFE'
scen = 's1'
# pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '160v100v80'
data = 'by_gene'
p_case = 160
p_case_fun = p_case_syn = p_int_fun = p_int_syn = p_exp = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 100
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
Nref = 500
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
pi_tar1 = 0.20 #pi.target for NFE: 0.20 for s1 or 0 for s2
pi_tar2 = 0.80 #pi.target for AFR: 0.80 for s1 or 1 for s2
# genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)


dir_leg = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/', folder, '/')
dir_in = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/', folder, '/datasets/', scen, '/')
dir_out = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Results/')
# dir_out = paste0('/home/math/siglersa/mastersProject/Output/', pruning, '/', data, '/')

dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')

# Vectors to store unadjusted p-values
prox_ext_genes_p = prox_int_genes_p = c() #proxECAT
prox_weighted_ext_genes_p = prox_weighted_int_genes_p = c() #proxECAT-weighted
prox2_ext_genes_p = prox2_int_genes_p = prox2_all_genes_p = c() #LogProx
iecat_genes_p = skato_ext_genes_p = skato_int_genes_p = skato_all_genes_p = c() #iECAT-O and SKAT-O
skat_ext_genes_p = skat_int_genes_p = skat_all_genes_p = c() #SKAT
burden_ext_genes_p = burden_int_genes_p = burden_all_genes_p = c() #Burden

# Vectors to store adjusted p-values
prox_ext_genes_p_adj = prox_int_genes_p_adj = c() #proxECAT
prox_weighted_ext_genes_p_adj = prox_weighted_int_genes_p_adj = c() #proxECAT-weighted
prox2_ext_genes_p_adj = prox2_int_genes_p_adj = prox2_all_genes_p_adj = c() #LogProx
iecat_genes_p_adj = skato_ext_genes_p_adj = skato_int_genes_p_adj = skato_all_genes_p_adj = c() #iECAT-O and SKAT-O
skat_ext_genes_p_adj = skat_int_genes_p_adj = skat_all_genes_p_adj = c() #SKAT
burden_ext_genes_p_adj = burden_int_genes_p_adj = burden_all_genes_p_adj = c() #Burden


# loop through the simulation replicates
set.seed(1) 
i=1
for (i in 1:100){
  
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
  hap_cases = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn)
  # hap_cases_power = read_hap_homo(dir_in, Pop, i, "cases", p_case_fun, p_case_syn) # pcase % fun 100% syn
  # hap_cases_t1e = read_hap_homo(dir_in, Pop, i, "cases", p_int_fun, p_int_syn) # 100% fun 100% syn
  hap_int = read_hap(dir_in, Pop1, Pop2, i, scen, "internal.controls", p_int_fun, p_int_syn)
  hap_cc = read_hap(dir_in, Pop1, Pop2, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  hap_ref_pop1 = read_ref(dir_in, Pop1, i, scen, p_exp, p_exp)
  hap_ref_pop2 = read_ref(dir_in, Pop2, i, scen, p_exp, p_exp)
  
  # FOR POWER ONLY
  # Create a new hap cases dataframe that merges the cases used for power and t1e but only contains
  # the genes associated with each calculation
  # hap_cases = merge_cases(hap_cases_power, hap_cases_t1e, leg, genes_power)
  
  # convert the haplotypes into genotypes
  geno_cases = make_geno(hap_cases)
  geno_int = make_geno(hap_int)
  geno_cc = make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_cases = calc_allele_freqs(geno_cases, Ncase)
  count_int = calc_allele_freqs(geno_int, Nint)
  count_cc = calc_allele_freqs(geno_cc, Ncc)
  count_all = calc_allele_freqs_all(count_cases, count_int, count_cc, Ncase, Nint, Ncc)
  
  count_ref_pop1 = calc_allele_freqs_ref(Pop1, hap_ref_pop1, Nref)
  count_ref_pop2 = calc_allele_freqs_ref(Pop2, hap_ref_pop2, Nref)
  
  cc_refs = cbind(count_cc, count_ref_pop1, count_ref_pop2)
  
  # Estimate ancestry proportions using only COMMON variants
  cc_est_prop = est_props(cc_refs, Pop1, Pop2, maf)
  # cases_est_prop = est_props(cases_refs, Pop1, Pop2, maf)
  # int_est_prop = est_props(int_refs, Pop1, Pop2, maf)
  # 
  # prop_ests_cc <- rbind(prop_ests_cc, cc_est_prop)
  # prop_ests_cases <- rbind(prop_ests_cases, cases_est_prop)
  # prop_ests_int <- rbind(prop_ests_int, int_est_prop)
  
  # Calculate adjusted AFs 
  count_cc_adj = calc_adjusted_AF(cc_refs, Pop1, Pop2, cc_est_prop, pi_tar1, pi_tar2, Nref, Ncc)
  
  # identify the common variants
  common_ext = leg[which(count_cases$maf > maf | count_cc$maf > maf),]
  common_all = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc$maf > maf),]
  common_int = leg[which(count_cases$maf > maf | count_int$maf > maf),]
  
  common_ext_adj = leg[which(count_cases$maf > maf | count_cc_adj$maf > maf),]
  common_all_adj = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc_adj$maf > maf),]
  
  # convert genotypes into long format for ProxECAT v2
  data_cases = make_long(count_cases, leg, "case", "int")
  data_int = make_long(count_int, leg, "control", "int")
  data_cc = make_long(count_cc, leg, "control", "ext")
  data_cc_adj = make_long_adj(count_cc_adj, leg, "control", "ext")
  
  # combine the data together
  data_all = data.frame(lapply(rbind(data_cases, data_int, data_cc), factor)) %>%
    filter(!(id %in% common_all$id))
  
  data_prox = data.frame(lapply(rbind(data_cases, data_cc), factor)) %>%
    filter(!(id %in% common_ext$id))
  
  data_int = data.frame(lapply(rbind(data_cases, data_int), factor)) %>%
    filter(!(id %in% common_int$id))
  
  data_all_adj = data.frame(lapply(rbind(data_cases, data_int, data_cc_adj), factor)) %>%
    filter(!(id %in% common_all_adj$id))
  
  data_prox_adj = data.frame(lapply(rbind(data_cases, data_cc_adj), factor)) %>%
    filter(!(id %in% common_ext_adj$id))
  
  # proxECAT
  counts_int_wide = prox_gene_data_prep(data_cases, data_int, common_int)
  counts_ext_wide = prox_gene_data_prep(data_cases, data_cc, common_ext)
  ounts_ext_wide_adj = prox_gene_data_prep(data_cases, data_cc_adj, common_ext_adj)
  
  # Store the proxECAT and proxECAT-weighted p-values
  prox_int_genes_p = rbind(prox_int_genes_p, counts_int_wide$prox)
  prox_ext_genes_p = rbind(prox_ext_genes_p, counts_ext_wide$prox)
  prox_ext_genes_p_adj = rbind(prox_ext_genes_p_adj, counts_ext_wide_adj$prox)
  
  prox_weighted_int_genes_p = rbind(prox_weighted_int_genes_p, counts_int_wide$prox_w)
  prox_weighted_ext_genes_p = rbind(prox_weighted_ext_genes_p, counts_ext_wide$prox_w)
  prox_weighted_ext_genes_p_adj = rbind(prox_weighted_ext_genes_p_adj, counts_ext_wide_adj$prox_w)
  
  # create case/control phenotype matrices for iECAT/SKAT
  pheno_int = rep(0, (ncol(geno_cases) + ncol(geno_int)))
  pheno_int[1:ncol(geno_cases)] = 1
  
  pheno_ext = rep(0, (ncol(geno_cases) + ncol(geno_cc)))
  pheno_ext[1:ncol(geno_cases)] = 1
  
  pheno_all = rep(0, (ncol(geno_cases) + ncol(geno_int) + ncol(geno_cc)))
  pheno_all[1:ncol(geno_cases)] = 1
  
  # null model object
  obj_int = SKAT_Null_Model(as.numeric(pheno_int) ~ 1, out_type="D") # D-dichotomous
  obj_ext = SKAT_Null_Model(as.numeric(pheno_ext) ~ 1, out_type="D") # D-dichotomous
  obj_all = SKAT_Null_Model(as.numeric(pheno_all) ~ 1, out_type="D") # D-dichotomous
  
  # create combined genotype matrices
  geno_cases_int = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_int$row),] #internal
  geno_cases_cc = cbind(geno_cases, geno_cc, gene=leg$gene)[-union(leg_syn$row, common_ext$row),] #external
  geno_all = cbind(geno_cases, geno_int, geno_cc, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #internal+external
  geno_int_all = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O

  geno_ext = cbind(count_cc, gene=leg$gene)[-union(leg_syn$row, common_all$row),]
  geno_ext_adj = cbind(count_cc_adj, gene=leg$gene)[-union(leg_syn$row, common_all_adj$row),]
  
  # some colnames are same from cbinding the geno matrices, need to make them unique
  colnames(geno_cases_int) <- make.unique(colnames(geno_cases_int))
  colnames(geno_cases_cc) <- make.unique(colnames(geno_cases_cc))
  colnames(geno_all) <- make.unique(colnames(geno_all))
  colnames(geno_int_all) <- make.unique(colnames(geno_int_all))
  
  # create MAC matrix for external controls
  tbl = data.frame(a0=geno_ext$mac) %>% mutate(a1=2*Ncc-a0, gene = geno_ext$gene)
  tbl_adj = data.frame(a0=geno_ext_adj$mac) %>% mutate(a1=2*Ncc-a0, gene = geno_ext_adj$gene)
  
  # call ProxECATv2/iECAT/SKAT once per gene
  prox2_int_genes = prox2_ext_genes = prox2_all_genes = c()
  iecat_genes = skato_int_genes = skato_ext_genes = skato_all_genes = c()
  skat_int_genes = skat_ext_genes = skat_all_genes = c()
  burden_int_genes = burden_ext_genes = burden_all_genes = c()
  
  prox2_ext_genes_adj = prox2_all_genes_adj = c()
  iecat_genes_adj = skato_ext_genes_adj = skato_all_genes_adj = c()
  skat_ext_genes_adj = skat_all_genes_adj = c()
  burden_ext_genes_adj = burden_all_genes_adj = c()
  
  genes = levels(droplevels(as.factor(leg$gene)))
  # g = 1
  # gene_counts = leg %>% count(gene)
  # loop through the genes
  for(g in 1:length(genes)){
    
    # print(paste0('current gene: ', genes[g], ' (', g, ' of ', length(genes), ')'))
    
    # LogProx
    prox2_int = logprox_gene_data_prep(data_int, genes[g], data.all=FALSE)
    prox2_ext = logprox_gene_data_prep(data_prox, genes[g], data.all=FALSE)
    prox2_all = logprox_gene_data_prep(data_all, genes[g], data.all=TRUE)
    
    prox2_ext_adj = logprox_gene_data_prep(data_prox_adj, genes[g], data.all=FALSE)
    prox2_all_adj = logprox_gene_data_prep(data_all_adj, genes[g], data.all=TRUE)
    
    # Save the LogProx p-values
    prox2_int_genes = c(prox2_int_genes, prox2_int)
    prox2_ext_genes = c(prox2_ext_genes, prox2_ext)
    prox2_all_genes = c(prox2_all_genes, prox2_all)
    
    prox2_ext_genes_adj = c(prox2_ext_genes_adj, prox2_ext_adj)
    prox2_all_genes_adj = c(prox2_all_genes_adj, prox2_all_adj)
    
    ### Prepare data for iECAT and SKAT methods
    Z_int = geno_cases_int %>% filter(gene == genes[g]) %>% select(-gene)
    Z_ext = geno_cases_cc %>% filter(gene == genes[g]) %>% select(-gene)
    Z_all = geno_all %>% filter(gene == genes[g]) %>% select(-gene)
    Z_int_all = geno_int_all %>% filter(gene == genes[g]) %>% select(-gene)
    
    # subset the MAC matrix for the external controls for iECAT
    tbl_gene = tbl %>% filter(gene == genes[g]) %>% select(-gene)
    tbl_gene_adj = tbl_adj %>% filter(gene == genes[g]) %>% select(-gene)
    
    # call the iECAT-O function
    re_gene = iECAT(t(Z_int_all), obj_int, as.matrix(tbl_gene), method="optimal")
    re_gene_adj = iECAT(t(Z_int_all), obj_int, as.matrix(tbl_gene_adj), method="optimal")
    
    # Save the iECAT-O p-values
    iecat_genes = c(iecat_genes, re_gene$p.value)
    iecat_genes_adj = c(iecat_genes_adj, re_gene_adj$p.value)
    
    # call the SKAT-O functions
    skato_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="SKATO") # SKAT-O external
    skato_int_gene = SKATBinary(t(Z_int), obj_int, method="SKATO") # SKAT-O internal
    skato_all_gene = SKATBinary(t(Z_all), obj_all, method="SKATO") # SKAT-O internal+external
    
    # Save SKAT-O p-values
    skato_ext_genes = c(skato_ext_genes, skato_ext_gene$p.value)
    skato_int_genes = c(skato_int_genes, skato_int_gene$p.value)
    skato_all_genes = c(skato_all_genes, skato_all_gene$p.value)
    
    # Call the SKAT functions
    skat_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="SKAT") # SKAT external
    skat_int_gene = SKATBinary(t(Z_int), obj_int, method="SKAT") # SKAT internal
    skat_all_gene = SKATBinary(t(Z_all), obj_all, method="SKAT") # SKAT external
    
    # Save the SKAT p-values
    skat_ext_genes = c(skat_ext_genes, skat_ext_gene$p.value)
    skat_int_genes = c(skat_int_genes, skat_int_gene$p.value)
    skat_all_genes = c(skat_all_genes, skat_all_gene$p.value)
    
    # Call the Burden functions
    burden_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="Burden") # Burden external
    burden_int_gene = SKATBinary(t(Z_int), obj_int, method="Burden") # Burden internal
    burden_all_gene = SKATBinary(t(Z_all), obj_all, method="Burden") # Burden external
    
    # Save the Burden p-values
    burden_ext_genes = c(burden_ext_genes, burden_ext_gene$p.value)
    burden_int_genes = c(burden_int_genes, burden_int_gene$p.value)
    burden_all_genes = c(burden_all_genes, burden_all_gene$p.value)
  }
  
  # store the LogProx, iECAT-O, SKAT gene p-values
  # Each col represents a gene and each row represents a sim rep
  prox2_int_genes_p = rbind(prox2_int_genes_p, prox2_int_genes)
  prox2_ext_genes_p = rbind(prox2_ext_genes_p, prox2_ext_genes)
  prox2_all_genes_p = rbind(prox2_all_genes_p, prox2_all_genes)
  
  prox2_ext_genes_p_adj = rbind(prox2_ext_genes_p_adj, prox2_ext_genes_adj)
  prox2_all_genes_p_adj = rbind(prox2_all_genes_p, prox2_all_genes_adj)
  
  iecat_genes_p = rbind(iecat_genes_p, iecat_genes)
  iecat_genes_p_adj = rbind(iecat_genes_p_adj, iecat_genes_adj)
  
  skato_int_genes_p = rbind(skato_int_genes_p, skato_int_genes)
  skato_ext_genes_p = rbind(skato_ext_genes_p, skato_ext_genes)
  skato_all_genes_p = rbind(skato_all_genes_p, skato_all_genes)
  
  skat_ext_genes_p = rbind(skat_ext_genes_p, skat_ext_genes)
  skat_int_genes_p = rbind(skat_int_genes_p, skat_int_genes)
  skat_all_genes_p = rbind(skat_all_genes_p, skat_all_genes)
  
  burden_ext_genes_p = rbind(burden_ext_genes_p, burden_ext_genes)
  burden_int_genes_p = rbind(burden_int_genes_p, burden_int_genes)
  burden_all_genes_p = rbind(burden_all_genes_p, burden_all_genes)
  
  print(i)
}

# Set col names to the genes
colnames(prox_ext_genes_p) = colnames(prox_int_genes_p) = genes
colnames(prox_weighted_ext_genes_p) = colnames(prox_weighted_int_genes_p) = genes
colnames(prox2_ext_genes_p) = colnames(prox2_all_genes_p) = colnames(prox2_int_genes_p) = genes
colnames(iecat_genes_p) = colnames(skato_int_genes_p) = colnames(skato_ext_genes_p) = colnames(skato_all_genes_p) = genes
colnames(skat_int_genes_p) = colnames(skat_ext_genes_p) = colnames(skat_all_genes_p) = genes
colnames(burden_int_genes_p) = colnames(burden_ext_genes_p) = colnames(burden_all_genes_p) = genes

# Set file path name
file_path = paste0(int_prune, "_v_", ext_prune, "_", Pop1, "_", Pop2, "_", scen, "_maf", maf, ".txt")

#ProxECAT
write.table(prox_ext_genes_p, paste0(dir_out, "T1e_gene_prox_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_int_genes_p, paste0(dir_out, "T1e_gene_prox_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_int_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_int_", file_path), quote=F, row.names=F, col.names=T)
# LogProx
write.table(prox2_ext_genes_p, paste0(dir_out, "T1e_gene_prox2_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p, paste0(dir_out, "T1e_gene_prox2_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_int_genes_p, paste0(dir_out, "T1e_gene_prox2_int_", file_path), quote=F, row.names=F, col.names=T)
# iECAT-O and SKAT-O
write.table(iecat_genes_p, paste0(dir_out, "T1e_gene_iecat_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(skato_int_genes_p, paste0(dir_out, "T1e_gene_skato_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(skato_ext_genes_p, paste0(dir_out, "T1e_gene_skato_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(skato_all_genes_p, paste0(dir_out, "T1e_gene_skato_all_", file_path), quote=F, row.names=F, col.names=T)
# SKAT
write.table(skat_int_genes_p, paste0(dir_out, "T1e_gene_skat_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(skat_ext_genes_p, paste0(dir_out, "T1e_gene_skat_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(skat_all_genes_p, paste0(dir_out, "T1e_gene_skat_all_", file_path), quote=F, row.names=F, col.names=T)
# Burden
write.table(burden_int_genes_p, paste0(dir_out, "T1e_gene_burden_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(burden_ext_genes_p, paste0(dir_out, "T1e_gene_burden_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(burden_all_genes_p, paste0(dir_out, "T1e_gene_burden_all_", file_path), quote=F, row.names=F, col.names=T)
