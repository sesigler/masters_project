# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(devtools) #need for installing Summix
library(proxecat)
library(SKAT)
library(iECAT)

source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/read_in_funcs.R")
source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/general_data_manip.R")
source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/methods_funcs.R")
source("https://raw.githubusercontent.com/hendriau/Summix/main/R/adjAF.R")
# source("https://raw.githubusercontent.com/hendriau/Summix/main/R/summix.R")
source("https://raw.githubusercontent.com/sesigler/Summix/main/R/summix.R")

Pop_admx <- 'AFR_NFE'  
Pops <- c('AFR', 'NFE')
admx_props <- c(80, 20)
scen <- 'xSCENx'
sub_scen <- 'xSUBx'
p_case <- 160
p_case_fun <- p_case_syn <- p_int_fun <- p_int_syn <- 100
p_cc_fun <- p_cc_syn <- 80
Ncase <- 2000 
Nic <- 2000
Ncc <- 10000  
Nrefs <- c(2000, 2000)
maf <- 0.001 
genes_power <- c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)

end <- mysim
start <- end-99


dir_leg <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/pruned_haps/')
dir_in <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/datasets/')
dir_out <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/power/')

# dir_leg = dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/')

# Vectors to store p-values
prox_ext_genes_p <- prox_ext_genes_p_adj_Ncc <- prox_ext_genes_p_adj_Neff <- c() 
proxW_ext_genes_p <- proxW_ext_genes_p_adj_Ncc <- proxW_ext_genes_p_adj_Neff <- c() 
prox2_ext_genes_p <- prox2_ext_genes_p_adj_Ncc <- prox2_ext_genes_p_adj_Neff <- c()
prox2_all_genes_p <- prox2_all_genes_p_adj_Ncc <- prox2_all_genes_p_adj_Neff <- c()
iecat_genes_p <- iecat_genes_p_adj_Ncc <- iecat_genes_p_adj_Neff <- c() 
skato_int_genes_p <- skat_int_genes_p <- burden_int_genes_p <- c()

# Vectors to store proxECAT and LogProx AC info
prox_ext_ACs <- prox_ext_adj_Ncc_ACs <- prox_ext_adj_Neff_ACs <- c()
prox2_ext_ACs <- prox2_ext_adj_Ncc_ACs <- prox2_ext_adj_Neff_ACs <- c()
prox2_all_ACs <- prox2_all_adj_Ncc_ACs <- prox2_all_adj_Neff_ACs <- c()

# Vectors to store effective sample sizes
neff_vec <- c()


# loop through the simulation replicates
# i=1
for (i in start:end){
  
  # read in the legend file
  # leg = read_leg_homo(dir_leg, Pop, i)
  leg <- read.table(paste0(dir_leg, 'chr19.block37.', Pop_admx, '.sim', i, '.', p_case, 'fun.100syn.legend'), header=T, sep='\t') #RAREsim v2.1.1 pruning only
  leg$row <- 1:nrow(leg)
  
  # Need to mutate so counts get added up correctly for ZNF333
  leg <- leg %>% mutate(gene = ifelse(gene == "ZNF333;ZNF333(NM_001352243:exon9:UTR5)", "ZNF333", gene))
  
  # subset the synonymous variants from the legend file
  leg_syn <- leg %>% filter(fun=="syn")
  leg_fun <- leg %>% filter(fun=="fun")
  
  # read in the haplotype files
  hap_cases_power <- read_hap(dir_in, Pop_admx, i, scen, "cases", p_fun = p_case, p_syn = 100) # pcase % fun 100% syn
  hap_cases_t1e <- read_hap(dir_in, Pop_admx, i, scen, "cases", p_fun = 100, p_syn = 100) # 100% fun 100% syn
  hap_ic <- read_hap(dir_in, Pop_admx, i, scen, "internal.controls", p_fun = 100, p_syn = 100)
  hap_cc <- read_hap(dir_in, Pop_admx, i, scen, "common.controls", p_fun = p_cc_fun, p_syn = p_cc_syn)

  # Create empty list to store reference haplotytpes
  hap_refs <- setNames(vector("list", length(Pops)), paste0("hap_ref_pop", 1:length(Pops)))
  
  # Read in the reference haplotypes
  for (j in seq_along(Pops)) {
    hap_refs[[j]] <- read_hap(dir_in, Pops[j], i, scen, "refs", p_fun = 100, p_syn = 100)
  }
  
  # FOR POWER ONLY
  # Create a new hap cases dataframe that merges the cases used for power and t1e but only for
  # the genes associated with each calculation
  hap_case <- merge_cases(hap_cases_power, hap_cases_t1e, leg, genes_power)
  
  # convert the haplotypes into genotypes
  geno_case <- make_geno(hap_case)
  geno_ic <- make_geno(hap_ic)
  geno_cc <- make_geno(hap_cc)
  
  # calculate the allele counts/frequencies
  count_case <- calc_allele_freqs(geno_case, Ncase, Pop=NULL)
  count_ic <- calc_allele_freqs(geno_ic, Nic, Pop=NULL)
  count_cc <- calc_allele_freqs(geno_cc, Ncc, Pop=NULL)
  
  # Create empty list to store reference allele counts/frequencies
  count_refs <- list()
  
  # calculate the reference allele counts/frequencies
  for (j in seq_along(Pops)) {
    count_refs[[j]] <- calc_allele_freqs(hap_refs[[j]], Nrefs[j], Pop=Pops[j])
  }
  
  # Commbine data with references for Summix
  cc_refs <- cbind(count_cc, do.call(cbind, count_refs))
  case_refs <- cbind(count_case, do.call(cbind, count_refs))
  
  # Estimate ancestry proportions using only COMMON variants
  cc_est_prop <- est_props(cc_refs, Pops, maf)
  case_est_prop <- est_props(case_refs, Pops, maf)
  # int_est_prop = est_props(int_refs, Pop1, Pop2, maf)
  
  # prop_ests_cc <- rbind(prop_ests_cc, cc_est_prop)
  # prop_ests_cases <- rbind(prop_ests_cases, cases_est_prop)
  # prop_ests_int <- rbind(prop_ests_int, int_est_prop)
  
  # Calculate adjusted AFs
  count_cc_adj_Ncc <- calc_adjusted_AF(cc_refs, Pops, case_est_prop, cc_est_prop, Nref=Nrefs, Ncc, Neff=FALSE)
  adj_Neff <- calc_adjusted_AF(cc_refs, Pops, case_est_prop, cc_est_prop, Nref=Nrefs, Ncc, Neff=TRUE)
  
  # return counts and effective sample size from Neff adjusted data
  count_cc_adj_Neff <- adj_Neff[[1]]
  Neff <- adj_Neff[[2]]
  
  # Save effective sample size
  neff_vec <- c(neff_vec, Neff)
  
  # Identify variants where AF >= 1-maf
  flip_int <- which(count_case$af >= 1-maf | count_ic$af >= 1-maf)
  flip_ext <- which(count_case$af >= 1-maf | count_cc$af >= 1-maf)
  flip_all <- which(count_case$af >= 1-maf | count_ic$af >= 1-maf | count_cc$af > 1-maf)
  
  flip_ext_adj_Ncc <- which(count_case$af >= 1-maf | count_cc_adj_Ncc$af >= 1-maf)
  flip_all_adj_Ncc <- which(count_case$af >= 1-maf | count_ic$af >= 1-maf | count_cc_adj_Ncc$af >= 1-maf)
  
  flip_ext_adj_Neff <- which(count_case$af >= 1-maf | count_cc_adj_Neff$af >= 1-maf)
  flip_all_adj_Neff <- which(count_case$af >= 1-maf | count_ic$af >= 1-maf | count_cc_adj_Neff$af >= 1-maf)
  
  # Flip the data at the variants identified above for all the different combination of datasets
  # If no variants need to be flipped, return the unaltered datasets
  # Cases and internal controls
  int_data <- flip_data(leg, flip_int, geno_case, count_case, Ncase, cntrl="int", geno_ic, count_ic, Nic, 
                       geno.cc=NULL, count.cc=NULL, count.cc.adj=NULL, Ncc=NULL, adj=FALSE)
  
  leg_int <- int_data[[1]]
  geno_case_int <- int_data[[2]]
  geno_ic_int <- int_data[[3]]
  count_case_int <- int_data[[5]]
  count_ic_int <- int_data[[6]]
  
  # Cases and external controls
  ext_data <- flip_data(leg, flip_ext, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                       geno_cc, count_cc, count.cc.adj=NULL, Ncc, adj=FALSE)
  
  leg_ext <- ext_data[[1]]
  geno_case_ext <- ext_data[[2]]
  geno_cc_ext <- ext_data[[4]]
  count_case_ext <- ext_data[[5]]
  count_cc_ext <- ext_data[[7]]
  
  # Cases, internal controls, and external controls
  all_data <- flip_data(leg, flip_all, geno_case, count_case, Ncase, cntrl="all", geno_ic, count_ic, Nic, 
                       geno_cc, count_cc, count.cc.adj=NULL, Ncc, adj=FALSE)
  
  leg_all <- all_data[[1]]
  geno_case_all <- all_data[[2]]
  geno_ic_all <- all_data[[3]]
  geno_cc_all <- all_data[[4]]
  count_case_all <- all_data[[5]]
  count_ic_all <- all_data[[6]]
  count_cc_all <- all_data[[7]]
  
  # Cases and Ncc adjusted external controls
  ext_adj_Ncc_data <- flip_data(leg, flip_ext_adj_Ncc, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                               geno.cc=NULL, count.cc=NULL, count_cc_adj_Ncc, Ncc, adj=TRUE)
  
  leg_ext_adj_Ncc <- ext_adj_Ncc_data[[1]]
  geno_case_ext_adj_Ncc <- ext_adj_Ncc_data[[2]]
  count_case_ext_adj_Ncc <- ext_adj_Ncc_data[[5]]
  count_cc_ext_adj_Ncc <- ext_adj_Ncc_data[[8]]
  
  # Cases, internal controls, and Ncc adjusted external controls
  all_adj_Ncc_data <- flip_data(leg, flip_all_adj_Ncc, geno_case, count_case, Ncase, cntrl="all", geno_ic, count_ic, Nic, 
                               geno.cc=NULL, count.cc=NULL, count_cc_adj_Ncc, Ncc, adj=TRUE)
  
  leg_all_adj_Ncc <- all_adj_Ncc_data[[1]]
  geno_case_all_adj_Ncc <- all_adj_Ncc_data[[2]]
  geno_ic_all_adj_Ncc <- all_adj_Ncc_data[[3]]
  count_case_all_adj_Ncc <- all_adj_Ncc_data[[5]]
  count_ic_all_adj_Ncc <- all_adj_Ncc_data[[6]]
  count_cc_all_adj_Ncc <- all_adj_Ncc_data[[8]]
  
  # Cases and Neff adjusted external controls
  ext_adj_Neff_data <- flip_data(leg, flip_ext_adj_Neff, geno_case, count_case, Ncase, cntrl="ext", geno.ic=NULL, count.ic=NULL, Nic=NULL, 
                                geno.cc=NULL, count.cc=NULL, count_cc_adj_Neff, Ncc=Neff, adj=TRUE)
  
  leg_ext_adj_Neff <- ext_adj_Neff_data[[1]]
  geno_case_ext_adj_Neff <- ext_adj_Neff_data[[2]]
  count_case_ext_adj_Neff <- ext_adj_Neff_data[[5]]
  count_cc_ext_adj_Neff <- ext_adj_Neff_data[[8]]
  
  # Cases, internal controls, and Neff adjusted external controls
  all_adj_Neff_data <- flip_data(leg, flip_all_adj_Neff, geno_case, count_case, Ncase, cntrl="all", geno_ic, count_ic, Nic, 
                                geno.cc=NULL, count.cc=NULL, count_cc_adj_Neff, Ncc=Neff, adj=TRUE)
  
  leg_all_adj_Neff <- all_adj_Neff_data[[1]]
  geno_case_all_adj_Neff <- all_adj_Neff_data[[2]]
  geno_ic_all_adj_Neff <- all_adj_Neff_data[[3]]
  count_case_all_adj_Neff <- all_adj_Neff_data[[5]]
  count_ic_all_adj_Neff <- all_adj_Neff_data[[6]]
  count_cc_all_adj_Neff <- all_adj_Neff_data[[8]]
  
  # identify the common variants
  common_int <- leg[which(count_case_int$af > maf | count_ic_int$af > maf),]
  common_ext <- leg[which(count_case_ext$af > maf | count_cc_ext$af > maf),]
  common_all <- leg[which(count_case_all$af > maf | count_ic_all$af > maf | count_cc_all$af > maf),]
  
  common_ext_adj_Ncc <- leg[which(count_case_ext_adj_Ncc$af > maf | count_cc_ext_adj_Ncc$af > maf),]
  common_all_adj_Ncc <- leg[which(count_case_all_adj_Ncc$af > maf | count_ic_all_adj_Ncc$af > maf | count_cc_all_adj_Ncc$af > maf),]
  
  common_ext_adj_Neff <- leg[which(count_case_ext_adj_Neff$af > maf | count_cc_ext_adj_Neff$af > maf),]
  common_all_adj_Neff <- leg[which(count_case_all_adj_Neff$af > maf | count_ic_all_adj_Neff$af > maf | count_cc_all_adj_Neff$af > maf),]
  
  # proxECAT
  counts_ext_wide <- prox_gene_data_prep(count_case_ext, count_cc_ext, leg_ext, common_ext)
  counts_ext_wide_adj_Ncc <- prox_gene_data_prep(count_case_ext_adj_Ncc, count_cc_ext_adj_Ncc, leg_ext_adj_Ncc, common_ext_adj_Ncc)
  counts_ext_wide_adj_Neff <- prox_gene_data_prep(count_case_ext_adj_Neff, count_cc_ext_adj_Neff, leg_ext_adj_Neff, common_ext_adj_Neff)
  
  # Store the proxECAT and proxECAT-weighted p-values
  prox_ext_genes_p <- rbind(prox_ext_genes_p, counts_ext_wide$prox)
  prox_ext_genes_p_adj_Ncc <- rbind(prox_ext_genes_p_adj_Ncc, counts_ext_wide_adj_Ncc$prox)
  prox_ext_genes_p_adj_Neff <- rbind(prox_ext_genes_p_adj_Neff, counts_ext_wide_adj_Neff$prox)
  
  proxW_ext_genes_p <- rbind(proxW_ext_genes_p, counts_ext_wide$prox_w)
  proxW_ext_genes_p_adj_Ncc <- rbind(proxW_ext_genes_p_adj_Ncc, counts_ext_wide_adj_Ncc$prox_w)
  proxW_ext_genes_p_adj_Neff <- rbind(proxW_ext_genes_p_adj_Neff, counts_ext_wide_adj_Neff$prox_w)
  
  # Save the proxECAT AC info, add rep column, and move it to front of df
  prox_ext_ACs <- rbind(prox_ext_ACs, counts_ext_wide %>% mutate(rep = i) %>% relocate(rep))
  prox_ext_adj_Ncc_ACs <- rbind(prox_ext_adj_Ncc_ACs, counts_ext_wide_adj_Ncc %>% mutate(rep = i) %>% relocate(rep))
  prox_ext_adj_Neff_ACs <- rbind(prox_ext_adj_Neff_ACs, counts_ext_wide_adj_Neff %>% mutate(rep = i) %>% relocate(rep))
  
  
  ### Prep data for other methods
  # convert genotypes into long format for ProxECAT v2, combine datasets, and remove common variants
  data_prox <- format_logprox_data(leg_ext, count_case_ext, count_cc_ext, control_type="ext", count.control2=NULL, common_ext, data.all=FALSE)
  data_prox_adj_Ncc <- format_logprox_data(leg_ext_adj_Ncc, count_case_ext_adj_Ncc, count_cc_ext_adj_Ncc, control_type="ext", count.control2=NULL, common_ext_adj_Ncc, data.all=FALSE)
  data_prox_adj_Neff <- format_logprox_data(leg_ext_adj_Neff, count_case_ext_adj_Neff, count_cc_ext_adj_Neff, control_type="ext", count.control2=NULL, common_ext_adj_Neff, data.all=FALSE)
  
  data_all <- format_logprox_data(leg_all, count_case_all, count_ic_all, control_type="int", count.control2=count_cc_all, common_all, data.all=TRUE)
  data_all_adj_Ncc <- format_logprox_data(leg_all_adj_Ncc, count_case_all_adj_Ncc, count_ic_all_adj_Ncc, control_type="int", count.control2=count_cc_all_adj_Ncc, common_all_adj_Ncc, data.all=TRUE)
  data_all_adj_Neff <- format_logprox_data(leg_all_adj_Neff, count_case_all_adj_Neff, count_ic_all_adj_Neff, control_type="int", count.control2=count_cc_all_adj_Neff, common_all_adj_Neff, data.all=TRUE)
  
  # create case/control phenotype matrices for iECAT/SKAT
  pheno_int <- rep(0, (ncol(geno_case_int) + ncol(geno_ic_int)))
  pheno_int[1:ncol(geno_case_int)] <- 1
  
  # null model object
  obj_int <- SKAT_Null_Model(as.numeric(pheno_int) ~ 1, out_type="D") # D-dichotomous
  
  # create combined genotype matrices
  geno_iecat_int <- cbind(geno_case_all, geno_ic_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  geno_iecat_int_adj_Ncc <- cbind(geno_case_all_adj_Ncc, geno_ic_all_adj_Ncc, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Ncc$row),] #iECAT-O
  geno_iecat_int_adj_Neff <- cbind(geno_case_all_adj_Neff, geno_ic_all_adj_Neff, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Neff$row),] #iECAT-O
  
  geno_iecat_ext <- cbind(count_cc_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  geno_iecat_ext_adj_Ncc <- cbind(count_cc_all_adj_Ncc, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Ncc$row),] #iECAT-O
  geno_iecat_ext_adj_Neff <- cbind(count_cc_all_adj_Neff, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Neff$row),] #iECAT-O
  
  geno_skat_int <- cbind(geno_case_int, geno_ic_int, gene=leg$gene)[-union(leg_syn$row, common_int$row),] #internal
  
  # some colnames are same from cbinding the geno matrices, need to make them unique
  colnames(geno_iecat_int) <- make.unique(colnames(geno_iecat_int))
  colnames(geno_iecat_int_adj_Ncc) <- make.unique(colnames(geno_iecat_int_adj_Ncc))
  colnames(geno_iecat_int_adj_Neff) <- make.unique(colnames(geno_iecat_int_adj_Neff))
  colnames(geno_skat_int) <- make.unique(colnames(geno_skat_int))
  
  # create MAC matrix for external controls
  tbl <- data.frame(a0=geno_iecat_ext$ac) %>% mutate(a1=2*Ncc-a0, gene = geno_iecat_ext$gene)
  tbl_adj_Ncc <- data.frame(a0=geno_iecat_ext_adj_Ncc$ac) %>% mutate(a1=2*Ncc-a0, gene = geno_iecat_ext_adj_Ncc$gene)
  tbl_adj_Neff <- data.frame(a0=geno_iecat_ext_adj_Neff$ac) %>% mutate(a1=2*Neff-a0, gene = geno_iecat_ext_adj_Neff$gene)
  
  # call ProxECATv2/iECAT/SKAT once per gene
  prox2_ext_genes <- prox2_ext_genes_adj_Ncc <- prox2_ext_genes_adj_Neff <- c()
  prox2_all_genes <- prox2_all_genes_adj_Ncc <- prox2_all_genes_adj_Neff <- c()
  iecat_genes <- iecat_genes_adj_Ncc <- iecat_genes_adj_Neff <- c()
  skato_int_genes <- skat_int_genes <- burden_int_genes <- c()
  
  # Create vectors to store ACs for LogProx
  prox2_ext_genes_ACs <- prox2_ext_adj_Ncc_genes_ACs <- prox2_ext_adj_Neff_genes_ACs <- c()
  prox2_all_genes_ACs <- prox2_all_adj_Ncc_genes_ACs <- prox2_all_adj_Neff_genes_ACs <- c()
  
  genes <- levels(droplevels(as.factor(leg$gene)))
  # g = 1
  # gene_counts = leg %>% count(gene)
  # loop through the genes
  for(g in 1:length(genes)){
    
    # print(paste0('current gene: ', genes[g], ' (', g, ' of ', length(genes), ')'))
    
    # LogProx
    prox2_ext <- logprox_gene_data_prep(data_prox, genes[g], data.all=FALSE)
    prox2_ext_adj_Ncc <- logprox_gene_data_prep(data_prox_adj_Ncc, genes[g], data.all=FALSE)
    prox2_ext_adj_Neff <- logprox_gene_data_prep(data_prox_adj_Neff, genes[g], data.all=FALSE)
    
    prox2_all <- logprox_gene_data_prep(data_all, genes[g], data.all=TRUE)
    prox2_all_adj_Ncc <- logprox_gene_data_prep(data_all_adj_Ncc, genes[g], data.all=TRUE)
    prox2_all_adj_Neff <- logprox_gene_data_prep(data_all_adj_Neff, genes[g], data.all=TRUE)
    
    # Save the LogProx p-values
    prox2_ext_genes <- c(prox2_ext_genes, prox2_ext[[1]])
    prox2_ext_genes_adj_Ncc <- c(prox2_ext_genes_adj_Ncc, prox2_ext_adj_Ncc[[1]])
    prox2_ext_genes_adj_Neff <- c(prox2_ext_genes_adj_Neff, prox2_ext_adj_Neff[[1]])
    
    prox2_all_genes <- c(prox2_all_genes, prox2_all[[1]])
    prox2_all_genes_adj_Ncc <- c(prox2_all_genes_adj_Ncc, prox2_all_adj_Ncc[[1]])
    prox2_all_genes_adj_Neff <- c(prox2_all_genes_adj_Neff, prox2_all_adj_Neff[[1]])
    
    # Save the LogProx ACs (saving the simulation replicate number, the gene, and the ACs)
    prox2_ext_genes_ACs <- rbind(prox2_ext_genes_ACs, prox2_ext[[2]] %>% mutate(rep = i, gene = genes[g]) %>% relocate(rep, gene))
    prox2_ext_adj_Ncc_genes_ACs <- rbind(prox2_ext_adj_Ncc_genes_ACs, prox2_ext_adj_Ncc[[2]] %>% mutate(rep = i, gene = genes[g]) %>% relocate(rep, gene))
    prox2_ext_adj_Neff_genes_ACs <- rbind(prox2_ext_adj_Neff_genes_ACs, prox2_ext_adj_Neff[[2]] %>% mutate(rep = i, gene = genes[g]) %>% relocate(rep, gene))
    
    prox2_all_genes_ACs <- rbind(prox2_all_genes_ACs, prox2_all[[2]] %>% mutate(rep = i, gene = genes[g]) %>% relocate(rep, gene))
    prox2_all_adj_Ncc_genes_ACs <- rbind(prox2_all_adj_Ncc_genes_ACs, prox2_all_adj_Ncc[[2]] %>% mutate(rep = i, gene = genes[g]) %>% relocate(rep, gene))
    prox2_all_adj_Neff_genes_ACs <- rbind(prox2_all_adj_Neff_genes_ACs, prox2_all_adj_Neff[[2]] %>% mutate(rep = i, gene = genes[g]) %>% relocate(rep, gene))
    
    ### Prepare data for iECAT and SKAT methods
    Z_iecat <- geno_iecat_int %>% filter(gene == genes[g]) %>% select(-gene) #iECAT
    Z_iecat_adj_Ncc <- geno_iecat_int_adj_Ncc %>% filter(gene == genes[g]) %>% select(-gene) #iECAT
    Z_iecat_adj_Neff <- geno_iecat_int_adj_Neff %>% filter(gene == genes[g]) %>% select(-gene) #iECAT
    
    Z_int <- geno_skat_int %>% filter(gene == genes[g]) %>% select(-gene) #SKAT-O, SKAT, Burden
    
    # subset the MAC matrix for the external controls for iECAT
    tbl_gene <- tbl %>% filter(gene == genes[g]) %>% select(-gene)
    tbl_gene_adj_Ncc <- tbl_adj_Ncc %>% filter(gene == genes[g]) %>% select(-gene)
    tbl_gene_adj_Neff <- tbl_adj_Neff %>% filter(gene == genes[g]) %>% select(-gene)
    
    # call the iECAT-O function
    re_gene <- iECAT(t(Z_iecat), obj_int, as.matrix(tbl_gene), method="optimal")
    re_gene_adj_Ncc <- iECAT(t(Z_iecat_adj_Ncc), obj_int, as.matrix(tbl_gene_adj_Ncc), method="optimal")
    re_gene_adj_Neff <- iECAT(t(Z_iecat_adj_Neff), obj_int, as.matrix(tbl_gene_adj_Neff), method="optimal")
    
    # Save the iECAT-O p-values
    iecat_genes <- c(iecat_genes, re_gene$p.value)
    iecat_genes_adj_Ncc <- c(iecat_genes_adj_Ncc, re_gene_adj_Ncc$p.value)
    iecat_genes_adj_Neff <- c(iecat_genes_adj_Neff, re_gene_adj_Neff$p.value)
    
    # call the SKAT-O, SKAT, and Burden functions
    skato_int_gene <- SKATBinary(t(Z_int), obj_int, method="SKATO") # SKAT-O internal
    skat_int_gene <- SKATBinary(t(Z_int), obj_int, method="SKAT") # SKAT internal
    burden_int_gene <- SKATBinary(t(Z_int), obj_int, method="Burden") # Burden internal
    
    # Save SKAT-O, SKAT, and Burden p-values
    skato_int_genes <- c(skato_int_genes, skato_int_gene$p.value)
    skat_int_genes <- c(skat_int_genes, skat_int_gene$p.value)
    burden_int_genes <- c(burden_int_genes, burden_int_gene$p.value)
    
  }
  
  # Store the LogProx and iECAT-O gene p-values
  # Each col represents a gene and each row represents a sim rep
  
  # LogProx p-values
  prox2_ext_genes_p <- rbind(prox2_ext_genes_p, prox2_ext_genes)
  prox2_ext_genes_p_adj_Ncc <- rbind(prox2_ext_genes_p_adj_Ncc, prox2_ext_genes_adj_Ncc)
  prox2_ext_genes_p_adj_Neff <- rbind(prox2_ext_genes_p_adj_Neff, prox2_ext_genes_adj_Neff)
  
  prox2_all_genes_p <- rbind(prox2_all_genes_p, prox2_all_genes)
  prox2_all_genes_p_adj_Ncc <- rbind(prox2_all_genes_p_adj_Ncc, prox2_all_genes_adj_Ncc)
  prox2_all_genes_p_adj_Neff <- rbind(prox2_all_genes_p_adj_Neff, prox2_all_genes_adj_Neff)
  
  # LogProx AC info
  prox2_ext_ACs <- rbind(prox2_ext_ACs, prox2_ext_genes_ACs)
  prox2_ext_adj_Ncc_ACs <- rbind(prox2_ext_adj_Ncc_ACs, prox2_ext_adj_Ncc_genes_ACs)
  prox2_ext_adj_Neff_ACs <- rbind(prox2_ext_adj_Neff_ACs, prox2_ext_adj_Neff_genes_ACs)
  
  prox2_all_ACs <- rbind(prox2_all_ACs, prox2_all_genes_ACs)
  prox2_all_adj_Ncc_ACs <- rbind(prox2_all_adj_Ncc_ACs, prox2_all_adj_Ncc_genes_ACs)
  prox2_all_adj_Neff_ACs <- rbind(prox2_all_adj_Neff_ACs, prox2_all_adj_Neff_genes_ACs)
  
  # iECAT-O p-values
  iecat_genes_p <- rbind(iecat_genes_p, iecat_genes)
  iecat_genes_p_adj_Ncc <- rbind(iecat_genes_p_adj_Ncc, iecat_genes_adj_Ncc)
  iecat_genes_p_adj_Neff <- rbind(iecat_genes_p_adj_Neff, iecat_genes_adj_Neff)
  
  # SKAT-O, SKAT, Burden p-values
  skato_int_genes_p <- rbind(skato_int_genes_p, skato_int_genes)
  skat_int_genes_p <- rbind(skat_int_genes_p, skat_int_genes)
  burden_int_genes_p <- rbind(burden_int_genes_p, burden_int_genes)
  
  print(i)
}

# Set col names to the genes
colnames(prox_ext_genes_p) <- colnames(prox_ext_genes_p_adj_Ncc) <- colnames(prox_ext_genes_p_adj_Neff) <- genes
colnames(proxW_ext_genes_p) <- colnames(proxW_ext_genes_p_adj_Ncc) <- colnames(proxW_ext_genes_p_adj_Neff) <- genes
colnames(prox2_ext_genes_p) <- colnames(prox2_ext_genes_p_adj_Ncc) <- colnames(prox2_ext_genes_p_adj_Neff) <- genes
colnames(prox2_all_genes_p) <- colnames(prox2_all_genes_p_adj_Ncc) <- colnames(prox2_all_genes_p_adj_Neff) <- genes
colnames(iecat_genes_p) <- colnames(iecat_genes_p_adj_Ncc) <- colnames(iecat_genes_p_adj_Neff) <- genes
colnames(skato_int_genes_p) <- colnames(skat_int_genes_p) <- colnames(burden_int_genes_p) <- genes

# Set file path name
file_path <- paste0(scen, "_", sub_scen, "_maf", maf, ".txt")

# Save the proportion estimates
# write.table(data.frame(prop_ests_cc), paste0(dir_out, "T1e_cc_prop_ests_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), quote=F, row.names=F)

# Save the effective sample sizes
write.csv(neff_vec, paste0(dir_out, "neff_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)

# Save the AC info
write.csv(prox_ext_ACs, paste0(dir_out, "ACs_prox_ext_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(prox_ext_adj_Ncc_ACs, paste0(dir_out, "ACs_prox_ext_adj_Ncc_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(prox_ext_adj_Neff_ACs, paste0(dir_out, "ACs_prox_ext_adj_Neff_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)

write.csv(prox2_ext_ACs, paste0(dir_out, "ACs_prox2_ext_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(prox2_ext_adj_Ncc_ACs, paste0(dir_out, "ACs_prox2_ext_adj_Ncc_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(prox2_ext_adj_Neff_ACs, paste0(dir_out, "ACs_prox2_ext_adj_Neff_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)

write.csv(prox2_all_ACs, paste0(dir_out, "ACs_prox2_all_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(prox2_all_adj_Ncc_ACs, paste0(dir_out, "ACs_prox2_all_adj_Ncc_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(prox2_all_adj_Neff_ACs, paste0(dir_out, "ACs_prox2_all_adj_Neff_", scen, "_", sub_scen, "_maf", maf, ".csv"), quote=F, row.names=F)

# ProxECAT
write.table(prox_ext_genes_p, paste0(dir_out, "Power_gene_prox_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# ProxECAT-weighted
write.table(proxW_ext_genes_p, paste0(dir_out, "Power_gene_prox_weighted_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox_weighted_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox_weighted_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# LogProx
write.table(prox2_ext_genes_p, paste0(dir_out, "Power_gene_prox2_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox2_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox2_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p, paste0(dir_out, "Power_gene_prox2_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_prox2_all_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p_adj_Neff, paste0(dir_out, "Power_gene_prox2_all_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# iECAT-O
write.table(iecat_genes_p, paste0(dir_out, "Power_gene_iecat_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(iecat_genes_p_adj_Ncc, paste0(dir_out, "Power_gene_iecat_all_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(iecat_genes_p_adj_Neff, paste0(dir_out, "Power_gene_iecat_all_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# SKAT-O, SKAT, Burden
write.table(skato_int_genes_p, paste0(dir_out, "Power_gene_skato_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(skat_int_genes_p, paste0(dir_out, "Power_gene_skat_int_", file_path), quote=F, row.names=F, col.names=T)
write.table(burden_int_genes_p, paste0(dir_out, "Power_gene_burden_int_", file_path), quote=F, row.names=F, col.names=T)
