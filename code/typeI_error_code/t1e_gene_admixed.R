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

source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/read_in_funcs.R")
source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/general_data_manip.R")
source("https://raw.githubusercontent.com/sesigler/masters_project/main/code/functions/methods_funcs.R")
source("https://raw.githubusercontent.com/hendriau/Summix/main/R/adjAF.R")
source("https://raw.githubusercontent.com/hendriau/Summix/main/R/summix.R")

Pop_admx <- 'LTX'  
Pops <- c('IAM', 'NFE', 'EAS', 'AFR')
admx_props <- c(47, 44, 5, 4)
scen <- 's2'
sub_scen <- 'default'
p_case <- 160
p_case_fun <- p_case_syn <- p_int_fun <- p_int_syn <- 100
p_cc_fun <- p_cc_syn <- 80
Ncase <- 2000 
Nic <- 2000
Ncc <- 10000  
Nrefs <- c(2000, 2000, 2000, 2000)
maf <- 0.001 

end <- 100 # change back to mysim
start <- end-99


dir_leg <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/pruned_haps/')
dir_in <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/datasets/')
dir_out <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/t1e/')

# dir_leg <- dir_in <- paste0('C:/Users/sagee/Documents/HendricksLab/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/')
# dir_in <- paste0('C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/', sim_params, '/')
# dir_out <- 'C:/Users/sagee/Documents/HendricksLab/admixed/Sim_42k/'

# Vectors to store p-values
prox_ext_genes_p <- prox_ext_genes_p_adj_Ncc <- prox_ext_genes_p_adj_Neff <- c()
proxW_ext_genes_p <- proxW_ext_genes_p_adj_Ncc <- proxW_ext_genes_p_adj_Neff <- c() 
prox2_ext_genes_p <- prox2_ext_genes_p_adj_Ncc <- prox2_ext_genes_p_adj_Neff <- c()
prox2_all_genes_p <- prox2_all_genes_p_adj_Ncc <- prox2_all_genes_p_adj_Neff <- c()
iecat_genes_p <- iecat_genes_p_adj_Ncc <- iecat_genes_p_adj_Neff <- c() 
skato_int_genes_p <- skato_ext_genes_p <- skato_all_genes_p <- c() 
skat_int_genes_p <- skat_ext_genes_p <- skat_all_genes_p <- c() 
burden_int_genes_p <- burden_ext_genes_p <- burden_all_genes_p <- c() 

# Vectors to store proxECAT and LogProx AC info
prox_ext_ACs <- prox_ext_adj_Ncc_ACs <- prox_ext_adj_Neff_ACs <- c()
prox2_ext_ACs <- prox2_ext_adj_Ncc_ACs <- prox2_ext_adj_Neff_ACs <- c()
prox2_all_ACs <- prox2_all_adj_Ncc_ACs <- prox2_all_adj_Neff_ACs <- c()

# Vectors to store effective sample sizes
neff_vec <- c()

# loop through the simulation replicates
set.seed(1) 
# i=55
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
  hap_case <- read_hap(dir_in, Pop_admx, i, scen, "cases", p_case_fun, p_case_syn)
  hap_ic <- read_hap(dir_in, Pop_admx, i, scen, "internal.controls", p_int_fun, p_int_syn)
  hap_cc <- read_hap(dir_in, Pop_admx, i, scen, "common.controls", p_cc_fun, p_cc_syn)
  
  # Create empty list to store reference haplotytpes
  hap_refs <- setNames(vector("list", length(Pops)), paste0("hap_ref_pop", 1:length(Pops)))
  
  # Read in the reference haplotypes
  for (j in seq_along(Pops)) {
    hap_refs[[j]] <- read_hap(dir_in, Pops[j], i, scen, "refs", p_fun = 100, p_syn = 100)
  }
  
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

  pheno_ext <- rep(0, (ncol(geno_case_ext) + ncol(geno_cc_ext)))
  pheno_ext[1:ncol(geno_case_ext)] <- 1

  pheno_all <- rep(0, (ncol(geno_case_all) + ncol(geno_ic_all) + ncol(geno_cc_all)))
  pheno_all[1:ncol(geno_case_all)] <- 1

  # null model object
  obj_int <- SKAT_Null_Model(as.numeric(pheno_int) ~ 1, out_type="D") # D-dichotomous
  obj_ext <- SKAT_Null_Model(as.numeric(pheno_ext) ~ 1, out_type="D") # D-dichotomous
  obj_all <- SKAT_Null_Model(as.numeric(pheno_all) ~ 1, out_type="D") # D-dichotomous

  # create combined genotype matrices
  geno_iecat_int <- cbind(geno_case_all, geno_ic_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  geno_iecat_int_adj_Ncc <- cbind(geno_case_all_adj_Ncc, geno_ic_all_adj_Ncc, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Ncc$row),] #iECAT-O
  geno_iecat_int_adj_Neff <- cbind(geno_case_all_adj_Neff, geno_ic_all_adj_Neff, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Neff$row),] #iECAT-O

  geno_iecat_ext <- cbind(count_cc_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
  geno_iecat_ext_adj_Ncc <- cbind(count_cc_all_adj_Ncc, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Ncc$row),] #iECAT-O
  geno_iecat_ext_adj_Neff <- cbind(count_cc_all_adj_Neff, gene=leg$gene)[-union(leg_syn$row, common_all_adj_Neff$row),] #iECAT-O

  geno_skat_int <- cbind(geno_case_int, geno_ic_int, gene=leg$gene)[-union(leg_syn$row, common_int$row),] #internal
  geno_skat_ext <- cbind(geno_case_ext, geno_cc_ext, gene=leg$gene)[-union(leg_syn$row, common_ext$row),] #external
  geno_skat_all <- cbind(geno_case_all, geno_ic_all, geno_cc_all, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #internal+external

  # some colnames are same from cbinding the geno matrices, need to make them unique
  colnames(geno_iecat_int) <- make.unique(colnames(geno_iecat_int))
  colnames(geno_iecat_int_adj_Ncc) <- make.unique(colnames(geno_iecat_int_adj_Ncc))
  colnames(geno_iecat_int_adj_Neff) <- make.unique(colnames(geno_iecat_int_adj_Neff))
  
  colnames(geno_skat_int) <- make.unique(colnames(geno_skat_int))
  colnames(geno_skat_ext) <- make.unique(colnames(geno_skat_ext))
  colnames(geno_skat_all) <- make.unique(colnames(geno_skat_all))

  # create MAC matrix for external controls
  tbl <- data.frame(a0=geno_iecat_ext$ac) %>% mutate(a1=2*Ncc-a0, gene = geno_iecat_ext$gene)
  tbl_adj_Ncc <- data.frame(a0=geno_iecat_ext_adj_Ncc$ac) %>% mutate(a1=2*Ncc-a0, gene = geno_iecat_ext_adj_Ncc$gene)
  tbl_adj_Neff <- data.frame(a0=geno_iecat_ext_adj_Neff$ac) %>% mutate(a1=2*Neff-a0, gene = geno_iecat_ext_adj_Neff$gene)

  # call ProxECATv2/iECAT/SKAT once per gene
  prox2_ext_genes <- prox2_ext_genes_adj_Ncc <- prox2_ext_genes_adj_Neff <- c()
  prox2_all_genes <- prox2_all_genes_adj_Ncc <- prox2_all_genes_adj_Neff <- c()
  iecat_genes <- iecat_genes_adj_Ncc <- iecat_genes_adj_Neff <- c()
  skato_int_genes <- skato_ext_genes <- skato_all_genes <- c()
  skat_int_genes <- skat_ext_genes <- skat_all_genes <- c()
  burden_int_genes <- burden_ext_genes <- burden_all_genes <- c()
  
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

    Z_int <- geno_skat_int %>% filter(gene == genes[g]) %>% select(-gene)
    Z_ext <- geno_skat_ext %>% filter(gene == genes[g]) %>% select(-gene)
    Z_all <- geno_skat_all %>% filter(gene == genes[g]) %>% select(-gene)

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

    # call the SKAT-O functions
    skato_int_gene <- SKATBinary(t(Z_int), obj_int, method="SKATO") # SKAT-O internal
    skato_ext_gene <- SKATBinary(t(Z_ext), obj_ext, method="SKATO") # SKAT-O external
    skato_all_gene <- SKATBinary(t(Z_all), obj_all, method="SKATO") # SKAT-O internal+external

    # Save SKAT-O p-values
    skato_int_genes <- c(skato_int_genes, skato_int_gene$p.value)
    skato_ext_genes <- c(skato_ext_genes, skato_ext_gene$p.value)
    skato_all_genes <- c(skato_all_genes, skato_all_gene$p.value)

    # Call the SKAT functions
    skat_int_gene <- SKATBinary(t(Z_int), obj_int, method="SKAT") # SKAT internal
    skat_ext_gene <- SKATBinary(t(Z_ext), obj_ext, method="SKAT") # SKAT external
    skat_all_gene <- SKATBinary(t(Z_all), obj_all, method="SKAT") # SKAT internal+external

    # Save the SKAT p-values
    skat_int_genes <- c(skat_int_genes, skat_int_gene$p.value)
    skat_ext_genes <- c(skat_ext_genes, skat_ext_gene$p.value)
    skat_all_genes <- c(skat_all_genes, skat_all_gene$p.value)

    # Call the Burden functions
    burden_int_gene <- SKATBinary(t(Z_int), obj_int, method="Burden") # Burden internal
    burden_ext_gene <- SKATBinary(t(Z_ext), obj_ext, method="Burden") # Burden external
    burden_all_gene <- SKATBinary(t(Z_all), obj_all, method="Burden") # Burden internal+external

    # Save the Burden p-values
    burden_int_genes <- c(burden_int_genes, burden_int_gene$p.value)
    burden_ext_genes <- c(burden_ext_genes, burden_ext_gene$p.value)
    burden_all_genes <- c(burden_all_genes, burden_all_gene$p.value)
  }

  # Store the LogProx, iECAT-O, SKAT gene p-values
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

  # SKAT-O p-values
  skato_int_genes_p <- rbind(skato_int_genes_p, skato_int_genes)
  skato_ext_genes_p <- rbind(skato_ext_genes_p, skato_ext_genes)
  skato_all_genes_p <- rbind(skato_all_genes_p, skato_all_genes)

  # SKAT p-values
  skat_int_genes_p <- rbind(skat_int_genes_p, skat_int_genes)
  skat_ext_genes_p <- rbind(skat_ext_genes_p, skat_ext_genes)
  skat_all_genes_p <- rbind(skat_all_genes_p, skat_all_genes)

  # Burden p-values
  burden_int_genes_p <- rbind(burden_int_genes_p, burden_int_genes)
  burden_ext_genes_p <- rbind(burden_ext_genes_p, burden_ext_genes)
  burden_all_genes_p <- rbind(burden_all_genes_p, burden_all_genes)

  print(i)
}

# Set col names to the genes
colnames(prox_ext_genes_p) <- colnames(prox_ext_genes_p_adj_Ncc) <- colnames(prox_ext_genes_p_adj_Neff) <- genes
colnames(proxW_ext_genes_p) <- colnames(proxW_ext_genes_p_adj_Ncc) <- colnames(proxW_ext_genes_p_adj_Neff) <- genes
colnames(prox2_ext_genes_p) <- colnames(prox2_ext_genes_p_adj_Ncc) <- colnames(prox2_ext_genes_p_adj_Neff) <- genes
colnames(prox2_all_genes_p) <- colnames(prox2_all_genes_p_adj_Ncc) <- colnames(prox2_all_genes_p_adj_Neff) <- genes
colnames(iecat_genes_p) <- colnames(iecat_genes_p_adj_Ncc) <- colnames(iecat_genes_p_adj_Neff) <- genes
colnames(skato_int_genes_p) <- colnames(skato_ext_genes_p) <- colnames(skato_all_genes_p) <- genes
colnames(skat_int_genes_p) <- colnames(skat_ext_genes_p) <- colnames(skat_all_genes_p) <- genes
colnames(burden_int_genes_p) <- colnames(burden_ext_genes_p) <- colnames(burden_all_genes_p) <- genes

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

# Save the p-values
# ProxECAT
write.table(prox_ext_genes_p, paste0(dir_out, "T1e_gene_prox_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj_Ncc, paste0(dir_out, "T1e_gene_prox_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox_ext_genes_p_adj_Neff, paste0(dir_out, "T1e_gene_prox_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# ProxECAT-weighted
write.table(proxW_ext_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_genes_p_adj_Ncc, paste0(dir_out, "T1e_gene_prox_weighted_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(proxW_ext_genes_p_adj_Neff, paste0(dir_out, "T1e_gene_prox_weighted_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# LogProx
write.table(prox2_ext_genes_p, paste0(dir_out, "T1e_gene_prox2_ext_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p_adj_Ncc, paste0(dir_out, "T1e_gene_prox2_ext_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_ext_genes_p_adj_Neff, paste0(dir_out, "T1e_gene_prox2_ext_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p, paste0(dir_out, "T1e_gene_prox2_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p_adj_Ncc, paste0(dir_out, "T1e_gene_prox2_all_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p_adj_Neff, paste0(dir_out, "T1e_gene_prox2_all_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# iECAT-O
write.table(iecat_genes_p, paste0(dir_out, "T1e_gene_iecat_all_", file_path), quote=F, row.names=F, col.names=T)
write.table(iecat_genes_p_adj_Ncc, paste0(dir_out, "T1e_gene_iecat_all_adj_Ncc_", file_path), quote=F, row.names=F, col.names=T)
write.table(iecat_genes_p_adj_Neff, paste0(dir_out, "T1e_gene_iecat_all_adj_Neff_", file_path), quote=F, row.names=F, col.names=T)
# SKAT-O
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

##############################################################################
# Check adjustments between unadj AF and case AF
# plot(cc_refs$af, count_cases$af, xlim=c(0, 0.001), ylim = c(0, 0.001))
# abline(0, 1)
# 
# CCC_dat = CCC(targetAFs, adjustedAFs)$rho.c
# 
# # and adj AF and case AF
# plot(adj_AF$adjusted.AF$adjustedAF, count_cases$af, xlim=c(0, 0.001), ylim = c(0, 0.001))
# abline(0, 1)
# 
# # Combine the data
# data_refs <- cbind(count_cases, cc_refs)
# colnames(data_refs)[1:4] <- c("ac_case", "af_case", "ac_cc", "af_cc")
# data_refs$id = leg$id
# data_refs_rare = data_refs %>% filter(af_cc != 0)
# 
# case0 = leg[which((count_cases$af == 0)),]
# cc_unadj0 = leg[which((count_cc$af != 0)),]
# cc_adj0 = leg[which((cc_refs$adj_maf != 0)),]
# 
# case_cc_unadj0 = leg[which((count_cases$af != 0) | (count_cc$af != 0)),]
# case_cc_adj0 = leg[which((count_cases$af != 0) | (cc_refs$adj_maf != 0)),]
# 
# common_ext = leg[which((count_cases$af > maf & count_cases$af < 1-maf) | (count_cc$af > maf & count_cc$af < 1-maf)),]
# common_ext_adj = leg[which((count_cases$af > maf & count_cases$af < 1-maf) | counts_adj$af > maf),]
# 
# data_refs_rare_unadj = data_refs %>% filter(!(id %in% common_ext$id))
# data_refs_rare_adj = data_refs %>% filter(!(id %in% common_ext_adj$id))
# 
# #hist unadj af
# library(ggpubr)
# library(ggplot2)
# library(DescTools)
# 
# ccc_unadj <- CCC(data_refs_rare_unadj$af_case, data_refs_rare_unadj$af_cc)
# lab <- paste("CCC: ", round(ccc_unadj$rho.c[,1], digits = 3))
# plot(data_refs_rare_unadj$af_case, data_refs_rare_unadj$af_cc, xlab = "Case AF", ylab = "Unadjusted Control AF")
# title(main = "Case AF vs unadj External Control AF-Rare Variants")
# abline(0, 1)
# text(x = 0.0008, y = 0.0001, labels = lab)
# 
# ccc_adj <- CCC(data_refs_rare_adj$af_case, data_refs_rare_adj$adj_maf)
# lab <- paste("CCC: ", round(ccc_adj$rho.c[,1], digits = 3))
# plot(data_refs_rare_adj$af_case, data_refs_rare_adj$adj_maf, xlab = "Case AF", ylab = "adjusted Control AF")
# title(main = "Case AF vs adj External Control AF-Rare Variants")
# abline(0, 1)
# text(x = 0.0008, y = 0.0001, labels = lab)
# 
# ccc_adj_AC <- CCC(data_refs_rare_adj$ac_case, data_refs_rare_adj$adj_mac)
# lab <- paste("CCC: ", round(ccc_adj_AC$rho.c[,1], digits = 3))
# plot(data_refs_rare_adj$ac_case, data_refs_rare_adj$adj_mac, xlab = "Case AC", ylab = "adjusted Control AC")
# title(main = "Case AC vs adj External Control AC-Rare Variants")
# abline(0, 1)
# text(x = 8, y = 1, labels = lab)
# 
# g4 <- ggplot(data = data_refs_rare_adj, aes(x = ac_case, y = adj_mac)) +
#   stat_sum(aes(size = factor(after_stat(n))), geom = "point") +
#   labs(y='adjusted Control AC', x='Case AC', title = 'Case AC vs adjusted Controls AC-Rare Variants')
# g4
# 
# plot(data_refs$af_afr, data_refs$af_nfe)
# abline(0, 1)
# 
# ccc_ref_unadj <- CCC(data_refs_rare_unadj$af_afr, data_refs_rare_unadj$af_nfe)$rho.c
# plot(data_refs_rare_unadj$af_afr, data_refs_rare_unadj$af_nfe)
# abline(0, 1)
# 
# ccc_ref_adj <- CCC(data_refs_rare_adj$af_afr, data_refs_rare_adj$af_nfe, conf.level = 0.95)
# lab <- paste("CCC: ", round(ccc_ref_adj$rho.c[,1], digits = 3))
# plot(data_refs_rare_adj$af_afr, data_refs_rare_adj$af_nfe, xlab = "Case AF", ylab = "adjusted Control AF")
# title(main = "Case AF vs adj External Control AF-Rare Variants")
# abline(0, 1)
# text(x = 0.0009, y = 0.0001, labels = lab)
# 
# 
# # Histograms of other data when case AF is 0
# h_case0_unadjcc <- ggplot(data_refs %>% filter(af_case == 0), aes(x=af_cc)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case0_unadjcc
# 
# h_case0_adjcc <- ggplot(data_refs %>% filter(af_case == 0), aes(x=adj_maf)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case0_adjcc
# 
# h_case0_afr <- ggplot(data_refs %>% filter(af_case == 0), aes(x=af_afr)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case0_afr
# 
# h_case0_nfe <- ggplot(data_refs %>% filter(af_case == 0), aes(x=af_nfe)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case0_nfe
# 
# h_case0_all <- ggarrange(h_case0_unadjcc, h_case0_adjcc, h_case0_afr, h_case0_nfe, ncol=2, nrow=2)
# plot1 <- annotate_figure(h_case0_all, top = text_grob("Histogram of Unadjusted and Adjusted External Control and Reference AFs when Case AF is 0", size = 14))
# ggsave(file = paste0(dir_out, "caseAF_0_hist_cc_ref_AFs.jpg"), plot = plot1, height = 8, width = 16, units = 'in')
# 
# # Histograms of other data when case AC is 1
# h_case1_unadjcc <- ggplot(data_refs %>% filter(ac_case == 1), aes(x=ac_cc)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case1_unadjcc
# 
# h_case1_adjcc <- ggplot(data_refs %>% filter(ac_case == 1), aes(x=adj_mac)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case1_adjcc
# 
# h_case1_afr <- ggplot(data_refs %>% filter(ac_case == 1), aes(x=ac_afr)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case1_afr
# 
# h_case1_nfe <- ggplot(data_refs %>% filter(ac_case == 1), aes(x=ac_nfe)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# h_case1_nfe
# 
# h_case1_all <- ggarrange(h_case1_unadjcc, h_case1_adjcc, h_case1_afr, h_case1_nfe, ncol=2, nrow=2)
# h_case1 <- annotate_figure(h_case1_all, top = text_grob("Histogram of Unadjusted and Adjusted External Control and Reference ACs when Case AC is 1", size = 14))
# ggsave(file = paste0(dir_out, "caseAC_1_hist_cc_ref_ACs.jpg"), plot = h_case1, height = 8, width = 16, units = 'in')
# 
# # Histograms of other data when case AC is 2
# h_case5_unadjcc <- ggplot(data_refs %>% filter(ac_case == 5), aes(x=ac_cc)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# 
# h_case5_adjcc <- ggplot(data_refs %>% filter(ac_case == 5), aes(x=adj_mac)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# 
# h_case5_afr <- ggplot(data_refs %>% filter(ac_case == 5), aes(x=ac_afr)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# 
# h_case5_nfe <- ggplot(data_refs %>% filter(ac_case == 5), aes(x=ac_nfe)) +
#   geom_histogram(bins=10, color = "black", fill = "white")
# 
# h_case5_all <- ggarrange(h_case5_unadjcc, h_case5_adjcc, h_case5_afr, h_case5_nfe, ncol=2, nrow=2)
# h_case5 <- annotate_figure(h_case5_all, top = text_grob("Histogram of Unadjusted and Adjusted External Control and Reference ACs when Case AC is 5", size = 14))
# ggsave(file = paste0(dir_out, "caseAC_5_hist_cc_ref_ACs.jpg"), plot = h_case5, height = 8, width = 16, units = 'in')
# 
# h_cases <- ggplot(data_refs_rare, aes(x=af_case)) +
#   geom_histogram(bins=20, color = "black", fill = "white") +
#   xlim(-0.0001, 0.001)
# 
# h_cc <- ggplot(data_refs_rare, aes(x=af_cc)) +
#   geom_histogram(bins=20, color = "black", fill = "white") +
#   xlim(-0.0001, 0.001)
# 
# h_cc_adj <- ggplot(data_refs_rare, aes(x=adj_af)) +
#   geom_histogram(bins=20, color = "black", fill = "white") +
#   xlim(-0.0001, 0.001)
# 
# hists <- ggarrange(h_cases, h_cc, h_cc_adj, ncol=2, nrow=2)
##############################################################################

### Check dist
# data.cases = count.cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case", group="int") %>%
#   filter(mac!=0) %>% select(-count)
# data.cases2 = count_cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case", group="int") %>%
#   filter(mac!=0)
# 
# data.controls2 = count_cc %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="control", group="ext") %>%
#   filter(mac!=0)
# 
# names <- c("id", "gene", "fun", "case", "group")
# 
# data.prox2 = data.frame(rbind(data.cases2, data.controls2)) %>% mutate(across(all_of(names), as.factor)) %>%
#   filter(!(id %in% common_ext2$id))
# 
# # counts_gene = data.prox %>% group_by(gene, case, fun) %>% summarise(n = sum(mac))
# counts_gene = data.prox2 %>% group_by(gene, case, fun) %>% summarise(n = sum(mac))
# 
# counts_wide = tidyr::pivot_wider(counts_gene, names_from=c(case, fun), values_from=n,
#                                  values_fill=0, names_sep="_")
# 
# # Calculate the ratios
# counts_wide = counts_wide %>% mutate(case_ratio = case_fun/case_syn,
#                                      control_ratio = control_fun/control_syn)
# 
# # Calculate medians
# median_case_ratio = median(counts_wide$case_ratio)
# median_control_ratio = median(counts_wide$control_ratio)
# 
# # Calculate the weighted values and the p-values
# counts_wide = counts_wide %>% mutate(case_fun_w = case_fun / median_case_ratio,
#                                      control_fun_w = control_fun / median_control_ratio) %>%
#   mutate(prox = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA,
#                        proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
#          prox_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA,
#                          proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
# 
# library(ggplot2)
# library(ggpubr)
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/')
# case_af <- ggplot(data.prox %>% filter(case == "case"), aes(x=af)) +
#   geom_histogram(color="black", fill="white")
# case_af
# 
# case_ac <- ggplot(data.prox %>% filter(case == "case"), aes(x=ac)) +
#   geom_histogram(color="black", fill="white")
# case_ac
# 
# cc_ac <- ggplot(data.prox %>% filter(case == "control"), aes(x=ac)) +
#   geom_histogram(color="black", fill="white")
# cc_ac
# 
# plots <- ggarrange(case_af, case_ac, cc_ac, ncol=2, nrow=2)
# plots
# ggsave(file = paste0(dir_out, 'rare_caseAF_caseAC_adjccAC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
#        plot = plots, height = 8, width = 16, units = 'in')

##############################################################################
# adj_AF <- adjAF(data = cc_refs,
#                 reference = c("af_afr", "af_nfe"),
#                 observed = "af",
#                 pi.target = c(case_est_prop[, "af_afr"], case_est_prop[, "af_nfe"]), 
#                 pi.observed = c(cc_est_prop[, "af_afr"], cc_est_prop[, "af_nfe"]),
#                 adj_method = "average",
#                 N_reference = c(Nref, Nref),
#                 N_observed = Ncc,
#                 filter = TRUE) 
# 
# # Add adj AF to data frame
# cc_refs$adj_af <- adj_AF$adjusted.AF$adjustedAF
# 
# # Calculate the MINOR adjusted AF
# cc_refs$adj_maf <- ifelse(cc_refs$adj_af > .5, 1-cc_refs$adj_af, cc_refs$adj_af)
# 
# # Calculate the adjusted Minor AC using the effective samples size
# cc_refs$adj_mac <- round(cc_refs$adj_maf*(2*Ncc))
# 
# # Return just the adjusted data
# counts_adj <- cc_refs[, c("adj_mac", "adj_maf")]
# 
# # Rename columns so they are same as other data frames
# colnames(counts_adj) <- c("ac", "af")
# 
# 
# common_ext_adj2 = leg[which((count_cases$af > maf & count_cases$af < 1-maf) | counts_adj$af > maf),]
# 
# data.cases = count_cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case") %>% 
#   filter(ac!=0)
# 
# data.controls = count_cc %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="control") %>% 
#   filter(ac!=0)
# 
# data.controls.adj = counts_adj %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="control") %>% 
#   filter(ac!=0)
# 
# names <- c("id", "gene", "fun", "case")
# 
# data.prox = data.frame(rbind(data.cases, data.controls)) %>% mutate(across(all_of(names), as.factor)) %>% 
#   filter(!(id %in% common_ext$id))
# 
# data.prox.adj = data.frame(rbind(data.cases, data.controls.adj)) %>% mutate(across(all_of(names), as.factor)) %>% 
#   filter(!(id %in% common_ext_adj2$id))
# 
# # counts_gene = data.prox %>% group_by(gene, case, fun) %>% summarise(n = sum(mac))
# counts_gene = data.prox %>% group_by(gene, case, fun) %>% summarise(n = sum(ac))
# counts_gene_adj = data.prox.adj %>% group_by(gene, case, fun) %>% summarise(n = sum(ac))
# 
# counts_wide = tidyr::pivot_wider(counts_gene, names_from=c(case, fun), values_from=n,
#                                  values_fill=0, names_sep="_")
# 
# counts_wide_adj = tidyr::pivot_wider(counts_gene_adj, names_from=c(case, fun), values_from=n,
#                                      values_fill=0, names_sep="_")
# 
# # Calculate the ratios
# counts_wide = counts_wide %>% mutate(case_ratio = case_fun/case_syn,
#                                      control_ratio = control_fun/control_syn)
# 
# counts_wide_adj = counts_wide_adj %>% mutate(case_ratio = case_fun/case_syn,
#                                              control_ratio = control_fun/control_syn)
# 
# # Calculate medians
# median_case_ratio = median(counts_wide$case_ratio)
# median_control_ratio = median(counts_wide$control_ratio)
# 
# median_case_ratio_adj = median(counts_wide_adj$case_ratio)
# median_control_ratio_adj = median(counts_wide_adj$control_ratio)
# 
# # Calculate the weighted values and the p-values
# counts_wide = counts_wide %>% mutate(case_fun_w = case_fun / median_case_ratio,
#                                      control_fun_w = control_fun / median_control_ratio) %>% 
#   mutate(prox = ifelse((case_fun + control_fun < 5) | (case_syn + control_syn < 5), NA, 
#                        proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
#          prox_w = ifelse((case_fun_w + control_fun_w < 5) | (case_syn + control_syn < 5), NA, 
#                          proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
# 
# counts_wide_adj = counts_wide_adj %>% mutate(case_fun_w = case_fun / median_case_ratio_adj,
#                                              control_fun_w = control_fun / median_control_ratio_adj) %>% 
#   mutate(prox = ifelse((case_fun + control_fun < 5) | (case_syn + control_syn < 5), NA, 
#                        proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
#          prox_w = ifelse((case_fun_w + control_fun_w < 5) | (case_syn + control_syn < 5), NA, 
#                          proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
# 
# prox_ext_genes_p = rbind(prox_ext_genes_p, counts_wide$prox)
# prox_ext_genes_p_adj = rbind(prox_ext_genes_p_adj, counts_wide_adj$prox)
# prox_weighted_ext_genes_p = rbind(prox_weighted_ext_genes_p, counts_wide$prox_w)
# prox_weighted_ext_genes_p_adj = rbind(prox_weighted_ext_genes_p_adj, counts_wide_adj$prox_w)
##############################################################################
# Plot code
# Plot AFR vs NFE AF rare variants
# cc_dat2 <- CCC(count_all2$af_afr, count_all2$af_nfe, conf.level = 0.95)
# lab2 <- paste("CCC: ", round(cc_dat2$rho.c[,1], digits = 3), " (95% CI ",
#               round(cc_dat2$rho.c[,2], digits = 3), " - ",
#               round(cc_dat2$rho.c[,3], digits = 3), ")", sep = "")
# plot(count_all2$af_afr, count_all2$af_nfe, xlab = "AFR AF", ylab = "NFE AF")
# title(main = "AFR vs NFE AF-Rare Variants")
# abline(0, 1)
# text(x = 0.002, y = 0.0065, labels = lab2)
# 
# # Plot case AF vs Ref AF rare variants
# plot(count_all2$af_afr, count_all2$af_case, xlab = "AFR AF", ylab = "Case AF")
# title(main = "Case vs AFR AF-Rare Variants")
# 
# plot(count_all2$af_nfe, count_all2$af_case, xlab = "NFE AF", ylab = "Case AF")
# title(main = "Case vs NFE AF-Rare Variants")
# 
# library(ggplot2)
# library(ggpubr)
# dir_out = 'C:/Users/sagee/Documents/HendricksLab/admixed/'
# g1 <- ggplot(data = count_all2, aes(x = af_afr, y = af_case)) +
#   stat_sum(aes(size = factor(after_stat(n))), geom = "point") +
#   labs(y='Case AF', x='AFR AF', title = 'Case vs AFR AF-Rare Variants')
# g1
# 
# g2 <- ggplot(data = count_all2, aes(x = af_nfe, y = af_case)) +
#   stat_sum(aes(size = factor(after_stat(n))), geom = "point") +
#   labs(y='Case AF', x='NFE AF', title = 'Case vs NFE AF-Rare Variants')
# g2
# 
# case_v_ref_AF <- ggarrange(g1, g2, ncol=2, nrow=1)
# ggsave(file = paste0(dir_out, 'case_v_ref_AF_rare_variants.jpg'),
#        plot = case_v_ref_AF, height = 8, width = 15, units = 'in')
# 
# g3 <- ggplot(data = count_all2, aes(x = af_afr, y = af_cc)) +
#   stat_sum(aes(size = factor(after_stat(n))), geom = "point") +
#   labs(y='Control AF', x='AFR AF', title = 'Control vs AFR AF-Rare Variants')
# g3
# 
# g4 <- ggplot(data = count_all2, aes(x = af_nfe, y = af_cc)) +
#   stat_sum(aes(size = factor(after_stat(n))), geom = "point") +
#   labs(y='Control AF', x='NFE AF', title = 'Control vs NFE AF-Rare Variants')
# g4
# control_v_ref_AF <- ggarrange(g3, g4, ncol=2, nrow=1)
# ggsave(file = paste0(dir_out, 'control_v_ref_AF_rare_variants.jpg'),
#        plot = control_v_ref_AF, height = 8, width = 15, units = 'in')
# 
# case_v_ref_AF <- ggarrange(g3, g4, ncol=2, nrow=1)
