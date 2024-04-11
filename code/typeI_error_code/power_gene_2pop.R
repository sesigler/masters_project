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
p_case_fun = 160
int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 80
Ncase = Nic = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
Nref_pop1 = 704
Nref_pop2 = 642
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
sim_params = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_', Pop1, 'ref', Nref_pop1, '_', Pop2, 'ref', Nref_pop2)
genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)

dir_leg = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/pruned_haps/')
dir_in = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', folder, '/', sim_params, '/', scen, '/')
dir_out = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Results/Sim_', Nsim, '/', sim_params, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir_out = paste0('/home/math/siglersa/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Results/Sim_', Nsim, '/', sim_params, '/prox_gene_adj_', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir_out = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Results/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/', sim_params, '/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/', sim_params, '/')
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/'

# Vectors to store unadjusted p-values
prox_int_genes_p = prox_ext_genes_p = c() #proxECAT
prox_weighted_int_genes_p = prox_weighted_ext_genes_p = c() #proxECAT-weighted
prox2_int_genes_p = prox2_ext_genes_p = prox2_all_genes_p = c() #LogProx
iecat_genes_p = c() #iECAT-O and SKAT-O

# Vectors to store adjusted p-values
prox_ext_genes_p_adj_Ncc = prox_ext_genes_p_adj_Neff = c() #proxECAT
prox_weighted_ext_genes_p_adj_Ncc = prox_weighted_ext_genes_p_adj_Neff = c() #proxECAT-weighted
prox2_ext_genes_p_adj_Ncc = prox2_all_genes_p_adj_Ncc = prox2_ext_genes_p_adj_Neff = prox2_all_genes_p_adj_Neff = c() #LogProx
iecat_genes_p_adj_Ncc = iecat_genes_p_adj_Neff = c() #iECAT-O

# Vector to save proportion estimates and macs/mafs
# prop_ests_cc = c()
# allele_data_unadj = allele_data_adj = allele_data_adj_unrounded = c()
# ac_af_data = c()


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
  hap_cases_power = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_case, p_syn = 100) # pcase % fun 100% syn
  hap_cases_t1e = read_hap(dir_in, Pop1, Pop2, i, scen, "cases", p_fun = 100, p_syn = 100) # 100% fun 100% syn
  hap_ic = read_hap(dir_in, Pop1, Pop2, i, scen, "internal.controls", p_fun = 100, p_syn = 100)
  hap_cc = read_hap(dir_in, Pop1, Pop2, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  hap_ref_pop1 = read_ref(dir_in, Pop1, i, scen, p_fun = 100, p_syn = 100)
  hap_ref_pop2 = read_ref(dir_in, Pop2, i, scen, p_fun = 100, p_syn = 100)
  
  # FOR POWER ONLY
  # Create a new hap cases dataframe that merges the cases used for power and t1e but only contains
  # the genes associated with each calculation
  hap_case = merge_cases(hap_cases_power, hap_cases_t1e, leg, genes_power)
  
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
  # int_est_prop = est_props(int_refs, Pop1, Pop2, maf)
  
  # prop_ests_cc <- rbind(prop_ests_cc, cc_est_prop)
  # prop_ests_cases <- rbind(prop_ests_cases, cases_est_prop)
  # prop_ests_int <- rbind(prop_ests_int, int_est_prop)
  
  # Add row index column to cc_refs in case summix removes variants during adjustment 
  cc_refs$row <- 1:nrow(cc_refs)
  
  # Calculate adjusted AFs
  count_cc_adj_Ncc = calc_adjusted_AF(cc_refs, Pop1, Pop2, case_est_prop, cc_est_prop, Nref=c(Nref_pop1, Nref_pop2), Ncc, Neff=FALSE)
  adj_Neff = calc_adjusted_AF(cc_refs, Pop1, Pop2, case_est_prop, cc_est_prop, Nref=c(Nref_pop1, Nref_pop2), Ncc, Neff=TRUE)
  
  # return counts and effective sample size from Neff adjusted data
  count_cc_adj_Neff = adj_Neff[[1]]
  Neff = adj_Neff[[2]]
  
  # Identify variants where AF >= 1-maf
  flip_int = which(count_case$af >= 1-maf | count_ic$af >= 1-maf)
  flip_ext = which(count_case$af >= 1-maf | count_cc$af >= 1-maf)
  flip_all = which(count_case$af >= 1-maf | count_ic$af >= 1-maf | count_cc$af > 1-maf)
  
  flip_ext_adj_Ncc = which(count_case$af >= 1-maf | count_cc_adj_Ncc$af >= 1-maf)
  flip_all_adj_Ncc = which(count_case$af >= 1-maf | count_ic$af >= 1-maf | count_cc_adj_Ncc$af >= 1-maf)
  
  flip_ext_adj_Neff = which(count_case$af >= 1-maf | count_cc_adj_Neff$af >= 1-maf)
  flip_all_adj_Neff = which(count_case$af >= 1-maf | count_ic$af >= 1-maf | count_cc_adj_Neff$af >= 1-maf)
  
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
  
  # Cases, internal controls, and external controls
  all_data = flip_data(leg, flip_all, geno_case, count_case, Ncase, cntrl="all", geno_ic, count_ic, Nic, 
                       geno_cc, count_cc, count.cc.adj=NULL, Ncc, adj=FALSE)
  
  leg_all = all_data[[1]]
  geno_case_all = all_data[[2]]
  geno_ic_all = all_data[[3]]
  geno_cc_all = all_data[[4]]
  count_case_all = all_data[[5]]
  count_ic_all = all_data[[6]]
  count_cc_all = all_data[[7]]
  
  # Cases and Ncc adjusted external controls
  ext_adj_Ncc_data = flip_data(leg, flip_ext_adj_Ncc, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                               geno.cc=NULL, count.cc=NULL, count_cc_adj_Ncc, Ncc, adj=TRUE)
  
  leg_ext_adj_Ncc = ext_adj_Ncc_data[[1]]
  geno_case_ext_adj_Ncc = ext_adj_Ncc_data[[2]]
  count_case_ext_adj_Ncc = ext_adj_Ncc_data[[5]]
  count_cc_ext_adj_Ncc = ext_adj_Ncc_data[[8]]
  
  # Cases, internal controls, and Ncc adjusted external controls
  all_adj_Ncc_data = flip_data(leg, flip_all_adj_Ncc, geno_case, count_case, Ncase, cntrl="all", geno_ic, count_ic, Nic, 
                               geno.cc=NULL, count.cc=NULL, count_cc_adj_Ncc, Ncc, adj=TRUE)
  
  leg_all_adj_Ncc = all_adj_Ncc_data[[1]]
  geno_case_all_adj_Ncc = all_adj_Ncc_data[[2]]
  geno_ic_all_adj_Ncc = all_adj_Ncc_data[[3]]
  count_case_all_adj_Ncc = all_adj_Ncc_data[[5]]
  count_ic_all_adj_Ncc = all_adj_Ncc_data[[6]]
  count_cc_all_adj_Ncc = all_adj_Ncc_data[[8]]
  
  # Cases and Neff adjusted external controls
  ext_adj_Neff_data = flip_data(leg, flip_ext_adj_Neff, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                                geno.cc=NULL, count.cc=NULL, count_cc_adj_Neff, Ncc=Neff, adj=TRUE)
  
  leg_ext_adj_Neff = ext_adj_Neff_data[[1]]
  geno_case_ext_adj_Neff = ext_adj_Neff_data[[2]]
  count_case_ext_adj_Neff = ext_adj_Neff_data[[5]]
  count_cc_ext_adj_Neff = ext_adj_Neff_data[[8]]
  
  # Cases, internal controls, and Neff adjusted external controls
  all_adj_Neff_data = flip_data(leg, flip_all_adj_Neff, geno_case, count_case, Ncase, cntrl="all", geno_ic, count_ic, Nic, 
                                geno.cc=NULL, count.cc=NULL, count_cc_adj_Neff, Ncc=Neff, adj=TRUE)
  
  leg_all_adj_Neff = all_adj_Neff_data[[1]]
  geno_case_all_adj_Neff = all_adj_Neff_data[[2]]
  geno_ic_all_adj_Neff = all_adj_Neff_data[[3]]
  count_case_all_adj_Neff = all_adj_Neff_data[[5]]
  count_ic_all_adj_Neff = all_adj_Neff_data[[6]]
  count_cc_all_adj_Neff = all_adj_Neff_data[[8]]
  
  # identify the common variants
  common_int = leg[which(count_case_int$af > maf | count_ic_int$af > maf),]
  common_ext = leg[which(count_case_ext$af > maf | count_cc_ext$af > maf),]
  common_all = leg[which(count_case_all$af > maf | count_ic_all$af > maf | count_cc_all$af > maf),]
  
  common_ext_adj_Ncc = leg[which(count_case_ext_adj_Ncc$af > maf | count_cc_ext_adj_Ncc$af > maf),]
  common_all_adj_Ncc = leg[which(count_case_all_adj_Ncc$af > maf | count_ic_all_adj_Ncc$af > maf | count_cc_all_adj_Ncc$af > maf),]
  
  common_ext_adj_Neff = leg[which(count_case_ext_adj_Neff$af > maf | count_cc_ext_adj_Neff$af > maf),]
  common_all_adj_Neff = leg[which(count_case_all_adj_Neff$af > maf | count_ic_all_adj_Neff$af > maf | count_cc_all_adj_Neff$af > maf),]
  
  # proxECAT
  counts_int_wide = prox_gene_data_prep(count_case_int, count_ic_int, leg_int, common_int)
  counts_ext_wide = prox_gene_data_prep(count_case_ext, count_cc_ext, leg_ext, common_ext)
  counts_ext_wide_adj_Ncc = prox_gene_data_prep(count_case_ext_adj_Ncc, count_cc_ext_adj_Ncc, leg_ext_adj_Ncc, common_ext_adj_Ncc)
  counts_ext_wide_adj_Neff = prox_gene_data_prep(count_case_ext_adj_Neff, count_cc_ext_adj_Neff, leg_ext_adj_Neff, common_ext_adj_Neff)
  
  # Store the proxECAT and proxECAT-weighted p-values
  prox_int_genes_p = rbind(prox_int_genes_p, counts_int_wide$prox)
  prox_ext_genes_p = rbind(prox_ext_genes_p, counts_ext_wide$prox)
  prox_ext_genes_p_adj_Ncc = rbind(prox_ext_genes_p_adj_Ncc, counts_ext_wide_adj_Ncc$prox)
  prox_ext_genes_p_adj_Neff = rbind(prox_ext_genes_p_adj_Neff, counts_ext_wide_adj_Neff$prox)
  
  prox_weighted_int_genes_p = rbind(prox_weighted_int_genes_p, counts_int_wide$prox_w)
  prox_weighted_ext_genes_p = rbind(prox_weighted_ext_genes_p, counts_ext_wide$prox_w)
  prox_weighted_ext_genes_p_adj_Ncc = rbind(prox_weighted_ext_genes_p_adj_Ncc, counts_ext_wide_adj_Ncc$prox_w)
  prox_weighted_ext_genes_p_adj_Neff = rbind(prox_weighted_ext_genes_p_adj_Neff, counts_ext_wide_adj_Neff$prox_w)
  
  ### Prep data for other methods
  # convert genotypes into long format for ProxECAT v2, combine datasets, and remove common variants
  data_int = format_logprox_data(leg_int, count_case_int, count_ic_int, control_type="int", count.control2=NULL, common_int, data.all=FALSE)
  
  data_prox = format_logprox_data(leg_ext, count_case_ext, count_cc_ext, control_type="ext", count.control2=NULL, common_ext, data.all=FALSE)
  data_prox_adj_Ncc = format_logprox_data(leg_ext_adj_Ncc, count_case_ext_adj_Ncc, count_cc_ext_adj_Ncc, control_type="ext", count.control2=NULL, common_ext_adj_Ncc, data.all=FALSE)
  data_prox_adj_Neff = format_logprox_data(leg_ext_adj_Neff, count_case_ext_adj_Neff, count_cc_ext_adj_Neff, control_type="ext", count.control2=NULL, common_ext_adj_Neff, data.all=FALSE)
  
  data_all = format_logprox_data(leg_all, count_case_all, count_ic_all, control_type="int", count.control2=count_cc_all, common_all, data.all=TRUE)
  data_all_adj_Ncc = format_logprox_data(leg_all_adj_Ncc, count_case_all_adj_Ncc, count_ic_all_adj_Ncc, control_type="int", count.control2=count_cc_all_adj_Ncc, common_all_adj_Ncc, data.all=TRUE)
  data_all_adj_Neff = format_logprox_data(leg_all_adj_Neff, count_case_all_adj_Neff, count_ic_all_adj_Neff, control_type="int", count.control2=count_cc_all_adj_Neff, common_all_adj_Neff, data.all=TRUE)
  
  # create case/control phenotype matrices for iECAT/SKAT
  pheno_int = rep(0, (ncol(geno_case_int) + ncol(geno_ic_int)))
  pheno_int[1:ncol(geno_case_int)] = 1
  
  # null model object
  obj_int = SKAT_Null_Model(as.numeric(pheno_int) ~ 1, out_type="D") # D-dichotomous
  
  # create combined genotype matrices
  geno_iecat_int = cbind(geno_case_all, geno_ic_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  geno_iecat_int_adj_Ncc = cbind(geno_case_all_adj_Ncc, geno_ic_all_adj_Ncc, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Ncc$row),] #iECAT-O
  geno_iecat_int_adj_Neff = cbind(geno_case_all_adj_Neff, geno_ic_all_adj_Neff, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Neff$row),] #iECAT-O
  
  geno_iecat_ext = cbind(count_cc_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  geno_iecat_ext_adj_Ncc = cbind(count_cc_all_adj_Ncc, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Ncc$row),] #iECAT-O
  geno_iecat_ext_adj_Neff = cbind(count_cc_all_adj_Neff, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Neff$row),] #iECAT-O
  
  # some colnames are same from cbinding the geno matrices, need to make them unique
  colnames(geno_iecat_int) <- make.unique(colnames(geno_iecat_int))
  colnames(geno_iecat_int_adj_Ncc) <- make.unique(colnames(geno_iecat_int_adj_Ncc))
  colnames(geno_iecat_int_adj_Neff) <- make.unique(colnames(geno_iecat_int_adj_Neff))
  
  # create MAC matrix for external controls
  tbl = data.frame(a0=geno_iecat_ext$ac) %>% mutate(a1=2*Ncc-a0, gene = geno_iecat_ext$gene)
  tbl_adj_Ncc = data.frame(a0=geno_iecat_ext_adj_Ncc$ac) %>% mutate(a1=2*Ncc-a0, gene = geno_iecat_ext_adj_Ncc$gene)
  tbl_adj_Neff = data.frame(a0=geno_iecat_ext_adj_Neff$ac) %>% mutate(a1=2*Neff-a0, gene = geno_iecat_ext_adj_Neff$gene)
  
  # call ProxECATv2/iECAT/SKAT once per gene
  prox2_int_genes = prox2_ext_genes = prox2_all_genes = c()
  iecat_genes = c()
  
  prox2_ext_genes_adj_Ncc = prox2_all_genes_adj_Ncc = prox2_ext_genes_adj_Neff = prox2_all_genes_adj_Neff = c()
  iecat_genes_adj_Ncc = iecat_genes_adj_Neff = c()
  
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
    
    prox2_ext_adj_Ncc = logprox_gene_data_prep(data_prox_adj_Ncc, genes[g], data.all=FALSE)
    prox2_all_adj_Ncc = logprox_gene_data_prep(data_all_adj_Ncc, genes[g], data.all=TRUE)
    
    prox2_ext_adj_Neff = logprox_gene_data_prep(data_prox_adj_Neff, genes[g], data.all=FALSE)
    prox2_all_adj_Neff = logprox_gene_data_prep(data_all_adj_Neff, genes[g], data.all=TRUE)
    
    # Save the LogProx p-values
    prox2_int_genes = c(prox2_int_genes, prox2_int)
    prox2_ext_genes = c(prox2_ext_genes, prox2_ext)
    prox2_all_genes = c(prox2_all_genes, prox2_all)
    
    prox2_ext_genes_adj_Ncc = c(prox2_ext_genes_adj_Ncc, prox2_ext_adj_Ncc)
    prox2_all_genes_adj_Ncc = c(prox2_all_genes_adj_Ncc, prox2_all_adj_Ncc)
    
    prox2_ext_genes_adj_Neff = c(prox2_ext_genes_adj_Neff, prox2_ext_adj_Neff)
    prox2_all_genes_adj_Neff = c(prox2_all_genes_adj_Neff, prox2_all_adj_Neff)
    
    
    ### Prepare data for iECAT and SKAT methods
    Z_iecat = geno_iecat_int %>% filter(gene == genes[g]) %>% select(-gene) #iECAT
    Z_iecat_adj_Ncc = geno_iecat_int_adj_Ncc %>% filter(gene == genes[g]) %>% select(-gene) #iECAT
    Z_iecat_adj_Neff = geno_iecat_int_adj_Neff %>% filter(gene == genes[g]) %>% select(-gene) #iECAT
    
    # subset the MAC matrix for the external controls for iECAT
    tbl_gene = tbl %>% filter(gene == genes[g]) %>% select(-gene)
    tbl_gene_adj_Ncc = tbl_adj_Ncc %>% filter(gene == genes[g]) %>% select(-gene)
    tbl_gene_adj_Neff = tbl_adj_Neff %>% filter(gene == genes[g]) %>% select(-gene)
    
    # call the iECAT-O function
    re_gene = iECAT(t(Z_iecat), obj_int, as.matrix(tbl_gene), method="optimal")
    re_gene_adj_Ncc = iECAT(t(Z_iecat_adj_Ncc), obj_int, as.matrix(tbl_gene_adj_Ncc), method="optimal")
    re_gene_adj_Neff = iECAT(t(Z_iecat_adj_Neff), obj_int, as.matrix(tbl_gene_adj_Neff), method="optimal")
    
    # Save the iECAT-O p-values
    iecat_genes = c(iecat_genes, re_gene$p.value)
    iecat_genes_adj_Ncc = c(iecat_genes_adj_Ncc, re_gene_adj_Ncc$p.value)
    iecat_genes_adj_Neff = c(iecat_genes_adj_Neff, re_gene_adj_Neff$p.value)
  }
  
  # store the LogProx, iECAT-O, SKAT gene p-values
  # Each col represents a gene and each row represents a sim rep
  prox2_int_genes_p = rbind(prox2_int_genes_p, prox2_int_genes)
  prox2_ext_genes_p = rbind(prox2_ext_genes_p, prox2_ext_genes)
  prox2_all_genes_p = rbind(prox2_all_genes_p, prox2_all_genes)
  
  prox2_ext_genes_p_adj_Ncc = rbind(prox2_ext_genes_p_adj_Ncc, prox2_ext_genes_adj_Ncc)
  prox2_all_genes_p_adj_Ncc = rbind(prox2_all_genes_p_adj_Ncc, prox2_all_genes_adj_Ncc)
  
  prox2_ext_genes_p_adj_Neff = rbind(prox2_ext_genes_p_adj_Neff, prox2_ext_genes_adj_Neff)
  prox2_all_genes_p_adj_Neff = rbind(prox2_all_genes_p_adj_Neff, prox2_all_genes_adj_Neff)
  
  iecat_genes_p = rbind(iecat_genes_p, iecat_genes)
  iecat_genes_p_adj_Ncc = rbind(iecat_genes_p_adj_Ncc, iecat_genes_adj_Ncc)
  iecat_genes_p_adj_Neff = rbind(iecat_genes_p_adj_Neff, iecat_genes_adj_Neff)
  
  print(i)
}

# Set col names to the genes
colnames(prox_int_genes_p) = colnames(prox_ext_genes_p) = genes
colnames(prox_weighted_int_genes_p) = colnames(prox_weighted_ext_genes_p) = genes
colnames(prox2_int_genes_p) = colnames(prox2_ext_genes_p) = colnames(prox2_all_genes_p) = genes
colnames(iecat_genes_p) = genes

colnames(prox_ext_genes_p_adj_Ncc) = colnames(prox_weighted_ext_genes_p_adj_Ncc) = genes
colnames(prox_ext_genes_p_adj_Neff) = colnames(prox_weighted_ext_genes_p_adj_Neff) = genes
colnames(prox2_ext_genes_p_adj_Ncc) = colnames(prox2_all_genes_p_adj_Ncc) = genes
colnames(prox2_ext_genes_p_adj_Neff) = colnames(prox2_all_genes_p_adj_Neff) = genes
colnames(iecat_genes_p_adj_Ncc) = colnames(iecat_genes_p_adj_Neff) = genes

# Set file path name
file_path = paste0(int_prune, "_v_", ext_prune, "_", Pop1, "_", Pop2, "_", scen, "_maf", maf, ".txt")

# Save the proportion estimates
# write.table(data.frame(prop_ests_cc), paste0(dir_out, "T1e_cc_prop_ests_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), quote=F, row.names=F)

# ProxECAT
write.table(prox_int_genes_p, paste0(dir_out, "Power_gene_prox_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p, paste0(dir_out, "Power_gene_prox_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# ProxECAT-weighted
write.table(prox_weighted_int_genes_p, paste0(dir_out, "Power_gene_prox_weighted_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p, paste0(dir_out, "Power_gene_prox_weighted_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox_weighted_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox_weighted_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# LogProx
write.table(prox2_int_genes_p, paste0(dir_out, "Power_gene_prox2_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p, paste0(dir_out, "Power_gene_prox2_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p, paste0(dir_out, "Power_gene_prox2_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox2_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox2_all_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox2_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox2_all_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# iECAT-O
write.table(iecat_genes_p, paste0(dir_out, "Power_gene_iecat_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(iecat_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_iecat_all_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(iecat_genes_p_adj_Neff, paste0(dir_out, "Power_gene_iecat_all_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)

