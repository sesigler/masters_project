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

source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/summix2_adjAF.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/summix2_summix.R")

# pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
# data = 'by_gene'
Pop1 = 'AFR'
Pop2 = 'NFE'
scen = 's2'
folder = '160v100v80'
p_case = 160
p_case_fun = p_case_syn = p_int_fun = p_int_syn = p_exp = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 100
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
Nref = 500
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
pi_tar1 = 1 #pi.target for AFR: 0.80 for s1 or 1 for s2
pi_tar2 = 0 #pi.target for NFE: 0.20 for s1 or 0 for s2
# genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)


dir_leg = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/', folder, '/')
dir_in = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/', folder, '/datasets/', scen, '/')
dir_out = paste0('/home/math/siglersa/admixed/', Pop1, '_', Pop2, '_pops/Results/')
# dir_out = paste0('/home/math/siglersa/mastersProject/Output/', pruning, '/', data, '/')

dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')

# Vectors to store unadjusted p-values
prox_int_genes_p = prox_ext_genes_p = c() #proxECAT
prox_weighted_int_genes_p = prox_weighted_ext_genes_p = c() #proxECAT-weighted
# prox2_int_genes_p = prox2_ext_genes_p = prox2_all_genes_p = c() #LogProx
# iecat_genes_p = c() #iECAT-O and SKAT-O
# skato_ext_genes_p = skato_int_genes_p = skato_all_genes_p = c() #iECAT-O and SKAT-O
# skat_ext_genes_p = skat_int_genes_p = skat_all_genes_p = c() #SKAT
# burden_ext_genes_p = burden_int_genes_p = burden_all_genes_p = c() #Burden

# Vectors to store adjusted p-values
prox_ext_genes_p_adj = c() #proxECAT
prox_weighted_ext_genes_p_adj = c() #proxECAT-weighted
# prox_ext_genes_p_adj_unrounded = prox_weighted_ext_genes_p_adj_unrounded = c() #proxECAT
# prox_ext_genes_p_adj_rounded = prox_weighted_ext_genes_p_adj_rounded = c() #proxECAT
# prox2_ext_genes_p_adj = prox2_all_genes_p_adj = c() #LogProx
# iecat_genes_p_adj = c() #iECAT-O and SKAT-O

# Vector to save proportion estimates and macs/mafs
# prop_ests_cc = c()
# allele_data_unadj = allele_data_adj = allele_data_adj_unrounded = c()
# ac_af_data = c()


# loop through the simulation replicates
set.seed(1) 
i=1
for (i in 1:10){
  
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
  # hap_cases_power = read_hap_homo(dir_in, Pop1, Pop2, i, scen, "cases", p_case_fun, p_case_syn) # pcase % fun 100% syn
  # hap_cases_t1e = read_hap_homo(dir_in, Pop1, Pop2, i, scen, "cases", p_int_fun, p_int_syn) # 100% fun 100% syn
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

  # prop_ests_cc <- rbind(prop_ests_cc, cc_est_prop)
  # prop_ests_cases <- rbind(prop_ests_cases, cases_est_prop)
  # prop_ests_int <- rbind(prop_ests_int, int_est_prop)
  
  ##############################################################################
  # Adjust AFs at the gene level instead of variant level
  # counts_cases_gene = count_cases %>% mutate(row = leg$row, gene = leg$gene, id = leg$id, fun = leg$fun)
  # cc_refs2 = cc_refs %>% mutate(row = leg$row, gene = leg$gene, id = leg$id, fun = leg$fun)
  # count_cc_gene = count_cc %>% mutate(row = leg$row, gene = leg$gene, id = leg$id, fun = leg$fun)
  # 
  # cc_refs_gene = cc_refs2 %>% filter(gene == "ZNF333")
  # counts_cases_gene = counts_cases_gene %>% filter(gene == "ZNF333")
  # count_cc_gene = count_cc_gene %>% filter(gene == "ZNF333")
  # 
  # # variants that are common in at least one dataset
  # common <- which(cc_refs_gene$maf > maf | cc_refs_gene[, "maf_afr"] > maf | cc_refs_gene[, "maf_nfe"] > maf) 
  # 
  # # Subset counts dataframe to only common variants
  # common_df <- cc_refs_gene[common,]
  # 
  # # Use summix to calculate ancestry proportion estimates
  # prop_est <- summix(data = common_df,
  #                    reference=c("maf_afr", #AFR
  #                                "maf_nfe"), #NFE
  #                    observed="maf") #leave out pi.start argument
  # 
  # adj_AF <- adjAF(data = cc_refs_gene,
  #                 reference = c("maf_afr", "maf_nfe"),
  #                 observed = "maf",
  #                 pi.target = c(pi_tar1, pi_tar2), #last one is AFR proportion
  #                 pi.observed = c(prop_est[, "maf_afr"], prop_est[, "maf_nfe"]),
  #                 adj_method = "average",
  #                 N_reference = c(Nref, Nref),
  #                 N_observed = Ncc,
  #                 filter = TRUE) 
  # 
  # # Add adj AF to data frame
  # cc_refs_gene$adj_maf <- adj_AF$adjusted.AF$adjustedAF
  # 
  # cc_refs_gene$adj_mac <- round(cc_refs_gene$adj_maf*(2*Ncc))
  # 
  # # Re-check that ACs and AFs are the minor ones (may be higher after adjustment)
  # counts_minor = calc_adj_allele_freqs(cc_refs_gene, Ncc)
  # 
  # # Create data frame with only the 2 adj MAC and AF columns
  # counts_adj <- counts_minor[, c("adj_mac2", "adj_maf2", "row", "gene", "id", "fun")]
  # 
  # # Rename columns so they are same as other data frames
  # colnames(counts_adj) <- c("mac", "maf", "row", "gene", "id", "fun")
  # 
  # common_ext = leg[which(counts_cases_gene$maf > maf | count_cc_gene$maf > maf),]
  # common_ext_adj = leg[which(counts_cases_gene$maf > maf | counts_adj$maf > maf),]
  # 
  # data.cases = counts_cases_gene %>% mutate(case="case", group="int") %>% 
  #   filter(mac!=0) %>% select(-count)
  # 
  # data.controls.adj = counts_adj %>% mutate(case="control", group="ext") %>% 
  #   filter(mac!=0)
  # 
  # data.controls = count_cc_gene %>% mutate(case="control", group="ext") %>% 
  #   filter(mac!=0) %>% select(-count)
  # 
  # names <- c("id", "gene", "fun", "case", "group")
  # 
  # data.prox.adj = data.frame(rbind(data.cases, data.controls.adj)) %>% mutate(across(all_of(names), as.factor)) %>% 
  #   filter(!(id %in% common_ext_adj$id))
  # 
  # data.prox = data.frame(rbind(data.cases, data.controls)) %>% mutate(across(all_of(names), as.factor)) %>% 
  #   filter(!(id %in% common_ext$id))
  # 
  # counts_gene_adj = data.prox.adj %>% group_by(gene, case, fun) %>% summarise(n = sum(mac))
  # 
  # counts_gene = data.prox %>% group_by(gene, case, fun) %>% summarise(n = sum(mac))
  # 
  # counts_wide_adj = tidyr::pivot_wider(counts_gene_adj, names_from=c(case, fun), values_from=n,
  #                                  values_fill=0, names_sep="_")
  # 
  # counts_wide = tidyr::pivot_wider(counts_gene, names_from=c(case, fun), values_from=n,
  #                                  values_fill=0, names_sep="_")
  # 
  # # Calculate the ratios
  # counts_wide_adj = counts_wide_adj %>% mutate(case_ratio = case_fun/case_syn,
  #                                              control_ratio = control_fun/control_syn)
  # 
  # counts_wide = counts_wide %>% mutate(case_ratio = case_fun/case_syn,
  #                                      control_ratio = control_fun/control_syn)
  # 
  # # Calculate medians
  # median_case_ratio_adj = median(counts_wide_adj$case_ratio)
  # median_control_ratio_adj = median(counts_wide_adj$control_ratio)
  # 
  # median_case_ratio = median(counts_wide$case_ratio)
  # median_control_ratio = median(counts_wide$control_ratio)
  # 
  # # Calculate the weighted values and the p-values
  # counts_wide_adj = counts_wide_adj %>% mutate(case_fun_w = case_fun / median_case_ratio_adj,
  #                                              control_fun_w = control_fun / median_control_ratio_adj) %>% 
  #   mutate(prox = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA, 
  #                        proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
  #          prox_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA, 
  #                          proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
  # 
  # counts_wide = counts_wide %>% mutate(case_fun_w = case_fun / median_case_ratio,
  #                                      control_fun_w = control_fun / median_control_ratio) %>% 
  #   mutate(prox = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA, 
  #                        proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
  #          prox_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA, 
  #                          proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
  ##############################################################################
  
  # Calculate adjusted AFs
  count_cc_adj = calc_adjusted_AF(cc_refs, Pop1, Pop2, cc_est_prop, pi_tar1, pi_tar2, Nref, Ncc)
  # count_cc_adj_unrounded = calc_adjusted_AF(cc_refs, Pop1, Pop2, cc_est_prop, pi_tar1, pi_tar2, Nref, Ncc, round_adj_mac = FALSE)

  # identify the common variants
  # common_int = leg[which(count_cases$af > maf & count_cases$af < 1-maf | count_int$af > maf & count_int$af < 1-maf),]
  common_ext = leg[which(count_cases$af > maf & count_cases$af < 1-maf | count_cc$af > maf & count_cc$af < 1-maf),]
  # common_all = leg[which(count_cases$af > maf & count_cases$af < 1-maf | count_int$af > maf & count_int$af < 1-maf | count_cc$af > maf & count_cc$af < 1-maf),]

  common_ext_adj = leg[which(count_cases$af > maf & count_cases$af < 1-maf | count_cc_adj$af > maf),]
  # common_all_adj = leg[which(count_cases$af > maf & count_cases$af < 1-maf | count_int$af > maf & count_int$af < 1-maf | count_cc_adj$af > maf),]
  
  ### Check dist
  # data.cases = count.cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case", group="int") %>% 
  #   filter(mac!=0) %>% select(-count)
  data.cases = count_cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case", group="int") %>% 
    filter(ac!=0)
  
  data.controls = count_cc_adj %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="control", group="ext") %>% 
    filter(ac!=0)
  
  names <- c("id", "gene", "fun", "case", "group")
  
  data.prox = data.frame(rbind(data.cases, data.controls)) %>% mutate(across(all_of(names), as.factor)) %>% 
    filter(!(id %in% common_ext_adj$id))
  
  library(ggplot2)
  library(ggpubr)
  dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
  case_af <- ggplot(data.prox %>% filter(case == "case"), aes(x=af)) +
    geom_histogram(color="black", fill="white")
  case_af
  
  case_ac <- ggplot(data.prox %>% filter(case == "case"), aes(x=ac)) +
    geom_histogram(color="black", fill="white")
  case_ac
  
  cc_ac <- ggplot(data.prox %>% filter(case == "control"), aes(x=ac)) +
    geom_histogram(color="black", fill="white")
  cc_ac
  
  plots <- ggarrange(case_af, case_ac, cc_ac, ncol=2, nrow=2)
  plots
  ggsave(file = paste0(dir_out, 'rare_caseAF_caseAC_adjccAC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
         plot = plots, height = 8, width = 16, units = 'in')

  # proxECAT
  # counts_int_wide = prox_gene_data_prep2(count_cases, count_int, leg, "int", adj=FALSE, common_int)
  counts_ext_wide = prox_gene_data_prep2(count_cases, count_cc, leg, "ext", adj=FALSE, common_ext)
  counts_ext_wide_adj = prox_gene_data_prep2(count_cases, count_cc_adj, leg, "ext", adj=TRUE, common_ext_adj)
  # counts_ext_wide_adj_unrounded = prox_gene_data_prep2(count_cases, count_cc_adj_unrounded, leg, "ext", adj=TRUE, common_ext_adj_unrounded, round_macs=FALSE)
  # counts_ext_wide_adj_rounded = prox_gene_data_prep2(count_cases, count_cc_adj_unrounded, leg, "ext", adj=TRUE, common_ext_adj_unrounded, round_macs=TRUE)

  # Store the proxECAT and proxECAT-weighted p-values
  # prox_int_genes_p = rbind(prox_int_genes_p, counts_int_wide$prox)
  prox_ext_genes_p = rbind(prox_ext_genes_p, counts_ext_wide$prox)
  prox_ext_genes_p_adj = rbind(prox_ext_genes_p_adj, counts_ext_wide_adj$prox)
  # prox_ext_genes_p_adj_unrounded = rbind(prox_ext_genes_p_adj_unrounded, counts_ext_wide_adj_unrounded$prox)
  # prox_ext_genes_p_adj_rounded = rbind(prox_ext_genes_p_adj_rounded, counts_ext_wide_adj_rounded$prox)

  # prox_weighted_int_genes_p = rbind(prox_weighted_int_genes_p, counts_int_wide$prox_w)
  prox_weighted_ext_genes_p = rbind(prox_weighted_ext_genes_p, counts_ext_wide$prox_w)
  prox_weighted_ext_genes_p_adj = rbind(prox_weighted_ext_genes_p_adj, counts_ext_wide_adj$prox_w)
  # prox_weighted_ext_genes_p_adj_unrounded = rbind(prox_weighted_ext_genes_p_adj_unrounded, counts_ext_wide_adj_unrounded$prox_w)
  # prox_weighted_ext_genes_p_adj_rounded = rbind(prox_weighted_ext_genes_p_adj_rounded, counts_ext_wide_adj_rounded$prox_w)

  # convert genotypes into long format for ProxECAT v2
  # data_cases = make_long(count_cases, leg, "case", "int")
  # data_int = make_long(count_int, leg, "control", "int")
  # data_cc = make_long(count_cc, leg, "control", "ext")
  # data_cc_adj = make_long_adj(count_cc_adj, leg, "control", "ext")
  # 
  # # combine the data together
  # data_int = data.frame(lapply(rbind(data_cases, data_int), factor)) %>%
  #   filter(!(id %in% common_int$id))
  # 
  # data_prox = data.frame(lapply(rbind(data_cases, data_cc), factor)) %>%
  #   filter(!(id %in% common_ext$id))
  # 
  # data_all = data.frame(lapply(rbind(data_cases, data_int, data_cc), factor)) %>%
  #   filter(!(id %in% common_all$id))
  # 
  # data_prox_adj = data.frame(lapply(rbind(data_cases, data_cc_adj), factor)) %>%
  #   filter(!(id %in% common_ext_adj$id))
  # 
  # data_all_adj = data.frame(lapply(rbind(data_cases, data_int, data_cc_adj), factor)) %>%
  #   filter(!(id %in% common_all_adj$id))
  # 
  # # create case/control phenotype matrices for iECAT/SKAT
  # pheno_int = rep(0, (ncol(geno_cases) + ncol(geno_int)))
  # pheno_int[1:ncol(geno_cases)] = 1
  #  
  # # pheno_ext = rep(0, (ncol(geno_cases) + ncol(geno_cc)))
  # # pheno_ext[1:ncol(geno_cases)] = 1
  # # 
  # # pheno_all = rep(0, (ncol(geno_cases) + ncol(geno_int) + ncol(geno_cc)))
  # # pheno_all[1:ncol(geno_cases)] = 1
  # 
  # # null model object
  # obj_int = SKAT_Null_Model(as.numeric(pheno_int) ~ 1, out_type="D") # D-dichotomous
  # # obj_ext = SKAT_Null_Model(as.numeric(pheno_ext) ~ 1, out_type="D") # D-dichotomous
  # # obj_all = SKAT_Null_Model(as.numeric(pheno_all) ~ 1, out_type="D") # D-dichotomous
  # 
  # # create combined genotype matrices
  # geno_int_all = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  # geno_int_all_adj = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_all_adj$row),] #iECAT-O
  # # geno_cases_int = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_int$row),] #internal
  # # geno_cases_cc = cbind(geno_cases, geno_cc, gene=leg$gene)[-union(leg_syn$row, common_ext$row),] #external
  # # geno_all = cbind(geno_cases, geno_int, geno_cc, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #internal+external
  # 
  # geno_ext = cbind(count_cc, gene=leg$gene)[-union(leg_syn$row, common_all$row),]
  # geno_ext_adj = cbind(count_cc_adj, gene=leg$gene)[-union(leg_syn$row, common_all_adj$row),]
  # 
  # # some colnames are same from cbinding the geno matrices, need to make them unique
  # colnames(geno_int_all) <- make.unique(colnames(geno_int_all))
  # colnames(geno_int_all_adj) <- make.unique(colnames(geno_int_all_adj))
  # # colnames(geno_cases_int) <- make.unique(colnames(geno_cases_int))
  # # colnames(geno_cases_cc) <- make.unique(colnames(geno_cases_cc))
  # # colnames(geno_all) <- make.unique(colnames(geno_all))
  # 
  # # create MAC matrix for external controls
  # tbl = data.frame(a0=geno_ext$mac) %>% mutate(a1=2*Ncc-a0, gene = geno_ext$gene)
  # tbl_adj = data.frame(a0=geno_ext_adj$mac) %>% mutate(a1=2*Ncc-a0, gene = geno_ext_adj$gene)
  # 
  # # call ProxECATv2/iECAT/SKAT once per gene
  # prox2_int_genes = prox2_ext_genes = prox2_all_genes = c()
  # iecat_genes = c()
  # # skato_int_genes = skato_ext_genes = skato_all_genes = c()
  # # skat_int_genes = skat_ext_genes = skat_all_genes = c()
  # # burden_int_genes = burden_ext_genes = burden_all_genes = c()
  # 
  # prox2_ext_genes_adj = prox2_all_genes_adj = c()
  # iecat_genes_adj = c()
  # 
  # genes = levels(droplevels(as.factor(leg$gene)))
  # # g = 1
  # # gene_counts = leg %>% count(gene)
  # # loop through the genes
  # for(g in 1:length(genes)){
  #   
  #   # print(paste0('current gene: ', genes[g], ' (', g, ' of ', length(genes), ')'))
  #   
  #   # LogProx
  #   prox2_int = logprox_gene_data_prep(data_int, genes[g], data.all=FALSE)
  #   prox2_ext = logprox_gene_data_prep(data_prox, genes[g], data.all=FALSE)
  #   prox2_all = logprox_gene_data_prep(data_all, genes[g], data.all=TRUE)
  #   
  #   prox2_ext_adj = logprox_gene_data_prep(data_prox_adj, genes[g], data.all=FALSE)
  #   prox2_all_adj = logprox_gene_data_prep(data_all_adj, genes[g], data.all=TRUE)
  #   
  #   # Save the LogProx p-values
  #   prox2_int_genes = c(prox2_int_genes, prox2_int)
  #   prox2_ext_genes = c(prox2_ext_genes, prox2_ext)
  #   prox2_all_genes = c(prox2_all_genes, prox2_all)
  #   
  #   prox2_ext_genes_adj = c(prox2_ext_genes_adj, prox2_ext_adj)
  #   prox2_all_genes_adj = c(prox2_all_genes_adj, prox2_all_adj)
  #   
  #   ### Prepare data for iECAT and SKAT methods
  #   Z_int_all = geno_int_all %>% filter(gene == genes[g]) %>% select(-gene)
  #   Z_int_all_adj = geno_int_all_adj %>% filter(gene == genes[g]) %>% select(-gene)
  #   # Z_int = geno_cases_int %>% filter(gene == genes[g]) %>% select(-gene)
  #   # Z_ext = geno_cases_cc %>% filter(gene == genes[g]) %>% select(-gene)
  #   # Z_all = geno_all %>% filter(gene == genes[g]) %>% select(-gene)
  #   
  #   # subset the MAC matrix for the external controls for iECAT
  #   tbl_gene = tbl %>% filter(gene == genes[g]) %>% select(-gene)
  #   tbl_gene_adj = tbl_adj %>% filter(gene == genes[g]) %>% select(-gene)
  #   
  #   # call the iECAT-O function
  #   re_gene = iECAT(t(Z_int_all), obj_int, as.matrix(tbl_gene), method="optimal")
  #   re_gene_adj = iECAT(t(Z_int_all_adj), obj_int, as.matrix(tbl_gene_adj), method="optimal")
  #   
  #   # Save the iECAT-O p-values
  #   iecat_genes = c(iecat_genes, re_gene$p.value)
  #   iecat_genes_adj = c(iecat_genes_adj, re_gene_adj$p.value)
  #   
  #   # # call the SKAT-O functions
  #   # skato_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="SKATO") # SKAT-O external
  #   # skato_int_gene = SKATBinary(t(Z_int), obj_int, method="SKATO") # SKAT-O internal
  #   # skato_all_gene = SKATBinary(t(Z_all), obj_all, method="SKATO") # SKAT-O internal+external
  #   # 
  #   # # Save SKAT-O p-values
  #   # skato_ext_genes = c(skato_ext_genes, skato_ext_gene$p.value)
  #   # skato_int_genes = c(skato_int_genes, skato_int_gene$p.value)
  #   # skato_all_genes = c(skato_all_genes, skato_all_gene$p.value)
  #   # 
  #   # # Call the SKAT functions
  #   # skat_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="SKAT") # SKAT external
  #   # skat_int_gene = SKATBinary(t(Z_int), obj_int, method="SKAT") # SKAT internal
  #   # skat_all_gene = SKATBinary(t(Z_all), obj_all, method="SKAT") # SKAT external
  #   # 
  #   # # Save the SKAT p-values
  #   # skat_ext_genes = c(skat_ext_genes, skat_ext_gene$p.value)
  #   # skat_int_genes = c(skat_int_genes, skat_int_gene$p.value)
  #   # skat_all_genes = c(skat_all_genes, skat_all_gene$p.value)
  #   # 
  #   # # Call the Burden functions
  #   # burden_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="Burden") # Burden external
  #   # burden_int_gene = SKATBinary(t(Z_int), obj_int, method="Burden") # Burden internal
  #   # burden_all_gene = SKATBinary(t(Z_all), obj_all, method="Burden") # Burden external
  #   # 
  #   # # Save the Burden p-values
  #   # burden_ext_genes = c(burden_ext_genes, burden_ext_gene$p.value)
  #   # burden_int_genes = c(burden_int_genes, burden_int_gene$p.value)
  #   # burden_all_genes = c(burden_all_genes, burden_all_gene$p.value)
  # }
  # 
  # # store the LogProx, iECAT-O, SKAT gene p-values
  # # Each col represents a gene and each row represents a sim rep
  # prox2_int_genes_p = rbind(prox2_int_genes_p, prox2_int_genes)
  # prox2_ext_genes_p = rbind(prox2_ext_genes_p, prox2_ext_genes)
  # prox2_all_genes_p = rbind(prox2_all_genes_p, prox2_all_genes)
  # 
  # prox2_ext_genes_p_adj = rbind(prox2_ext_genes_p_adj, prox2_ext_genes_adj)
  # prox2_all_genes_p_adj = rbind(prox2_all_genes_p, prox2_all_genes_adj)
  # 
  # iecat_genes_p = rbind(iecat_genes_p, iecat_genes)
  # iecat_genes_p_adj = rbind(iecat_genes_p_adj, iecat_genes_adj)
  
  # skato_int_genes_p = rbind(skato_int_genes_p, skato_int_genes)
  # skato_ext_genes_p = rbind(skato_ext_genes_p, skato_ext_genes)
  # skato_all_genes_p = rbind(skato_all_genes_p, skato_all_genes)
  # 
  # skat_ext_genes_p = rbind(skat_ext_genes_p, skat_ext_genes)
  # skat_int_genes_p = rbind(skat_int_genes_p, skat_int_genes)
  # skat_all_genes_p = rbind(skat_all_genes_p, skat_all_genes)
  # 
  # burden_ext_genes_p = rbind(burden_ext_genes_p, burden_ext_genes)
  # burden_int_genes_p = rbind(burden_int_genes_p, burden_int_genes)
  # burden_all_genes_p = rbind(burden_all_genes_p, burden_all_genes)
  
  print(i)
}

# Set col names to the genes
# colnames(prox_int_genes_p) = colnames(prox_ext_genes_p) = genes
colnames(prox_ext_genes_p) = genes
# colnames(prox_weighted_int_genes_p) = colnames(prox_weighted_ext_genes_p) = genes
colnames(prox_weighted_ext_genes_p) = genes
# colnames(prox2_int_genes_p) = colnames(prox2_ext_genes_p) = colnames(prox2_all_genes_p) = genes
# colnames(iecat_genes_p) = genes
# 
colnames(prox_ext_genes_p_adj) = colnames(prox_weighted_ext_genes_p_adj) = genes
# colnames(prox_ext_genes_p_adj_unrounded) = colnames(prox_weighted_ext_genes_p_adj_unrounded) = genes
# colnames(prox_ext_genes_p_adj_rounded) = colnames(prox_weighted_ext_genes_p_adj_rounded) = genes
# colnames(prox2_ext_genes_p_adj) = colnames(prox2_all_genes_p_adj) = genes
# colnames(iecat_genes_p_adj) = genes
# colnames(skato_int_genes_p) = colnames(skato_ext_genes_p) = colnames(skato_all_genes_p) = genes
# colnames(skat_int_genes_p) = colnames(skat_ext_genes_p) = colnames(skat_all_genes_p) = genes
# colnames(burden_int_genes_p) = colnames(burden_ext_genes_p) = colnames(burden_all_genes_p) = genes

# Set file path name
file_path = paste0(int_prune, "_v_", ext_prune, "_", Pop1, "_", Pop2, "_", scen, "_maf", maf, ".txt")

# colnames(allele_data) = c("case_mac", "case_maf", "unadj_cc_mac", "unadj_cc_maf", "adj_cc_mac", "adj_cc_maf", "gene")

# Save the proportion estimates
# write.table(data.frame(prop_ests_cc), paste0(dir_out, "T1e_cc_prop_ests_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), quote=F, row.names=F)
# write.table(data.frame(allele_data_unadj), paste0(dir_out, "T1e_MACs_MAFs_unadj_", file_path), quote=F, row.names=F)
# write.table(data.frame(allele_data_adj), paste0(dir_out, "T1e_MACs_MAFs_adj_", file_path), quote=F, row.names=F)
# write.table(data.frame(allele_data_adj_unrounded), paste0(dir_out, "T1e_MACs_MAFs_adj_unrounded_", file_path), quote=F, row.names=F)
# write.table(data.frame(ac_af_data), paste0(dir_out, "ac_af_data_", file_path), quote=F, row.names=F)

# ProxECAT
# write.table(prox_int_genes_p, paste0(dir_out, "T1e_gene_prox_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p, paste0(dir_out, "T1e_gene_prox_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj, paste0(dir_out, "T1e_gene_prox_ext_adj_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox_ext_genes_p_adj_unrounded, paste0(dir_out, "T1e_gene_prox_ext_adj_unrounded_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox_ext_genes_p_adj_rounded, paste0(dir_out, "T1e_gene_prox_ext_adj_rounded_", file_path), quote=F, row.names=F, col.names=T)
# ProxECAT-weighted
# write.table(prox_weighted_int_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p_adj, paste0(dir_out, "T1e_gene_prox_weighted_ext_adj_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox_weighted_ext_genes_p_adj_unrounded, paste0(dir_out, "T1e_gene_prox_weighted_ext_adj_unrounded_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox_weighted_ext_genes_p_adj_rounded, paste0(dir_out, "T1e_gene_prox_weighted_ext_adj_rounded_", file_path), quote=F, row.names=F, col.names=T)
# # LogProx
# write.table(prox2_int_genes_p, paste0(dir_out, "T1e_gene_prox2_int_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox2_ext_genes_p, paste0(dir_out, "T1e_gene_prox2_ext_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox2_all_genes_p, paste0(dir_out, "T1e_gene_prox2_all_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox2_ext_genes_p_adj, paste0(dir_out, "T1e_gene_prox2_ext_adj_", file_path), quote=F, row.names=F, col.names=T)
# write.table(prox2_all_genes_p_adj, paste0(dir_out, "T1e_gene_prox2_all_adj_", file_path), quote=F, row.names=F, col.names=T)
# # iECAT-O 
# write.table(iecat_genes_p, paste0(dir_out, "T1e_gene_iecat_all_", file_path), quote=F, row.names=F, col.names=T)
# write.table(iecat_genes_p_adj, paste0(dir_out, "T1e_gene_iecat_all_adj_", file_path), quote=F, row.names=F, col.names=T)

# # SKAT-O
# write.table(skato_int_genes_p, paste0(dir_out, "T1e_gene_skato_int_", file_path), quote=F, row.names=F, col.names=T)
# write.table(skato_ext_genes_p, paste0(dir_out, "T1e_gene_skato_ext_", file_path), quote=F, row.names=F, col.names=T)
# write.table(skato_all_genes_p, paste0(dir_out, "T1e_gene_skato_all_", file_path), quote=F, row.names=F, col.names=T)
# # SKAT
# write.table(skat_int_genes_p, paste0(dir_out, "T1e_gene_skat_int_", file_path), quote=F, row.names=F, col.names=T)
# write.table(skat_ext_genes_p, paste0(dir_out, "T1e_gene_skat_ext_", file_path), quote=F, row.names=F, col.names=T)
# write.table(skat_all_genes_p, paste0(dir_out, "T1e_gene_skat_all_", file_path), quote=F, row.names=F, col.names=T)
# # Burden
# write.table(burden_int_genes_p, paste0(dir_out, "T1e_gene_burden_int_", file_path), quote=F, row.names=F, col.names=T)
# write.table(burden_ext_genes_p, paste0(dir_out, "T1e_gene_burden_ext_", file_path), quote=F, row.names=F, col.names=T)
# write.table(burden_all_genes_p, paste0(dir_out, "T1e_gene_burden_all_", file_path), quote=F, row.names=F, col.names=T)
