# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(proxecat)
library(SKAT)
library(iECAT)

source("/home/math/siglersa/mastersProject/Input/read_in_funcs.R")
source("/home/math/siglersa/mastersProject/Input/general_data_manip.R")
source("/home/math/siglersa/mastersProject/Input/methods_funcs.R")

# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")

Pop = 'AFR'
# pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '160v100v80'
data = 'by_gene'
p_case = p_case_fun = 160
p_exp = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 80
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)
# scen = 's1'

# dir = '/storage/math/projects/compinfo/simulations/'
# setwd(paste0(dir, 'output/20K_', Pop, '/'))

# dir_leg = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/attempt2_combine_MACbins_legFiles_differ/') #only for 100v80
# dir_leg = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/')
# dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/')
# dir_out = paste0('/home/math/siglersa/mastersProject/Output/', pruning, '/', data, '/', folder, '/')
dir_leg = paste0('/home/math/siglersa/mastersProject/20K_', Pop, '/', folder, '/')
dir_in = paste0('/home/math/siglersa/mastersProject/20K_', Pop, '/', folder, '/')
dir_out = paste0('/home/math/siglersa/mastersProject/Output/20K_', Pop, '/', data, '/', folder, '/')
# dir_out = paste0('/home/math/siglersa/mastersProject/Output/', pruning, '/', data, '/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')

prox_ext_genes_p = prox_int_genes_p = c() #proxECAT
prox_weighted_ext_genes_p = prox_weighted_int_genes_p = c() #proxECAT-weighted
prox2_ext_genes_p = prox2_int_genes_p = prox2_all_genes_p = c() #LogProx
iecat_genes_p = skato_ext_genes_p = skato_int_genes_p = skato_all_genes_p = c() #iECAT-O and SKAT-O
skat_ext_genes_p = skat_int_genes_p = skat_all_genes_p = c() #SKAT
burden_ext_genes_p = burden_int_genes_p = burden_all_genes_p = c() #Burden


# loop through the simulation replicates
set.seed(1) 
# i=1
for (i in 1:100){
  
   # read in the legend file
   # leg = read_leg_homo(dir_leg, Pop, i)
   leg = read.table(paste0(dir_leg, 'chr19.block37.', Pop, '.sim', i, '.', p_case, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t') #RAREsim v2.1.1 pruning only
   leg$row = 1:nrow(leg)
   
   # Need to mutate so counts get added up correctly for ZNF333
   leg = leg %>% mutate(gene = ifelse(gene == "ZNF333;ZNF333(NM_001352243:exon9:UTR5)", "ZNF333", gene))
   
   # subset the synonymous variants from the legend file
   leg_syn = leg %>% filter(fun=="syn")
   leg_fun = leg %>% filter(fun=="fun")
   
   # read in the haplotype files
   # hap_cases = read_hap_homo(dir_in, Pop, i, "cases", p_case_fun, p_case_syn)
   hap_cases_power = read_hap_homo(dir_in, Pop, i, "cases", p_case_fun, p_case_syn) # pcase % fun 100% syn
   hap_cases_t1e = read_hap_homo(dir_in, Pop, i, "cases", p_int_fun, p_int_syn) # 100% fun 100% syn
   hap_int = read_hap_homo(dir_in, Pop, i, "internal.controls", p_int_fun, p_int_syn)
   hap_cc = read_hap_homo(dir_in, Pop, i, "common.controls", p_cc_fun, p_cc_syn)
   
   # Create a new hap cases dataframe that merges the cases used for power and t1e but only contains
   # the genes associated with each calculation
   hap_cases = merge_cases(hap_cases_power, hap_cases_t1e, leg, genes_power)

   # convert the haplotypes into genotypes
   geno_cases = make_geno(hap_cases)
   geno_int = make_geno(hap_int)
   geno_cc = make_geno(hap_cc)
   
   # calculate the allele counts/frequencies
   count_cases = calc_allele_freqs(geno_cases, Ncase)
   count_int = calc_allele_freqs(geno_int, Nint)
   count_cc = calc_allele_freqs(geno_cc, Ncc)
   count_all = calc_allele_freqs_all(count_cases, count_int, count_cc, Ncase, Nint, Ncc)
   
   # identify the common variants
   common_ext = leg[which(count_cases$maf > maf | count_cc$maf > maf),]
   common_all = leg[which(count_cases$maf > maf | count_int$maf > maf | count_cc$maf > maf),]
   common_int = leg[which(count_cases$maf > maf | count_int$maf > maf),]
   
   # convert genotypes into long format for ProxECAT v2
   data_cases = make_long(count_cases, leg, "case", "int")
   data_int = make_long(count_int, leg, "control", "int")
   data_cc = make_long(count_cc, leg, "control", "ext")

   # combine the data together
   data_all = data.frame(lapply(rbind(data_cases, data_int, data_cc), factor)) %>%
     filter(!(id %in% common_all$id))
   
   data_prox = data.frame(lapply(rbind(data_cases, data_cc), factor)) %>%
     filter(!(id %in% common_ext$id))
   
   data_int = data.frame(lapply(rbind(data_cases, data_int), factor)) %>%
     filter(!(id %in% common_int$id))
   
   # count the number of alleles per gene per status (case/control & fun/syn)
   # counts_gene = data_all %>% filter(!(case=="control" & group=="int")) %>% count(gene, case, fun)
   
   # proxECAT External
   counts_ext_gene = data_prox %>% count(gene, case, fun)
   counts_ext_wide = tidyr::pivot_wider(counts_ext_gene, names_from=c(case, fun), values_from=n,
                                        values_fill=0, names_sep="_")
   counts_ext_wide = counts_ext_wide %>% mutate(case_ratio = case_fun/case_syn,
                                                control_ratio = control_fun/control_syn,
                                                case_fun_w = case_fun/median(case_ratio),
                                                control_fun_w = control_fun/median(control_ratio)) %>%
     mutate(prox_ext = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA,
                              proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
            prox_ext_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA,
                                proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))

   # proxECAT Internal
   counts_int_gene = data_int %>% count(gene, case, fun)
   counts_int_wide = tidyr::pivot_wider(counts_int_gene, names_from=c(case, fun), values_from=n,
                                        values_fill=0, names_sep="_")
   counts_int_wide = counts_int_wide %>% mutate(case_ratio = case_fun/case_syn,
                                                control_ratio = control_fun/control_syn,
                                                case_fun_w = case_fun/median(case_ratio),
                                                control_fun_w = control_fun/median(control_ratio)) %>%
     mutate(prox_int = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA,
                              proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
            prox_int_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA,
                                proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))

   # Store the proxECAT and proxECAT-weighted p-values
   prox_ext_genes_p = rbind(prox_ext_genes_p, counts_ext_wide$prox_ext)
   prox_int_genes_p = rbind(prox_int_genes_p, counts_int_wide$prox_int)
   prox_weighted_ext_genes_p = rbind(prox_weighted_ext_genes_p, counts_ext_wide$prox_ext_w)
   prox_weighted_int_genes_p = rbind(prox_weighted_int_genes_p, counts_int_wide$prox_int_w)
     
   # counts_all = colSums(counts_wide[,-1])
   # names(counts_all) = colnames(counts_wide)[-1]

   # counts_all = data_all %>% count(case, fun)
   # counts_prox = data_prox %>% count(case, fun)
   # counts_int = data_int %>% count(case, fun)
   
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
   # geno_int_all = cbind(geno_cases, geno_int)[-union(leg_syn$row, common_all$row),] # for iECAT
   # geno_cases_int = cbind(geno_cases, geno_int)[-union(leg_syn$row, common_int$row),] # SKAT int
   # geno_cases_cc = cbind(geno_cases, geno_cc)[-union(leg_syn$row, common_ext$row),] # SKAT ext
   # geno_all = cbind(geno_cases, geno_int, geno_cc)[-union(leg_syn$row, common_all$row),] # SKAT (all)
   geno_cases_int = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_int$row),] #internal
   geno_cases_cc = cbind(geno_cases, geno_cc, gene=leg$gene)[-union(leg_syn$row, common_ext$row),] #external
   geno_all = cbind(geno_cases, geno_int, geno_cc, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #internal+external
   geno_int_all = cbind(geno_cases, geno_int, gene=leg$gene)[-union(leg_syn$row, common_all$row),] #iECAT-O
   # geno_ext = geno_cc[-union(leg_syn$row, common_all$row),] #will need to change geno_cc to counts_cc for admixed
   # geno_ext = count_cc[-union(leg_syn$row, common_all$row),] # need MAC for tbl object for iECAT
   geno_ext = cbind(count_cc, gene=leg$gene)[-union(leg_syn$row, common_all$row),]
   
   # some colnames are same from cbinding the geno matrices, need to make them unique
   colnames(geno_cases_int) <- make.unique(colnames(geno_cases_int))
   colnames(geno_cases_cc) <- make.unique(colnames(geno_cases_cc))
   colnames(geno_all) <- make.unique(colnames(geno_all))
   colnames(geno_int_all) <- make.unique(colnames(geno_int_all))

   # subset the genotype matrices
   # LAST FEW ROWS ARE NA FOR THE FIRST GENE, NEED TO FIGURE OUT WHY
   # I think it's because leg_fun still has the common variants in it
   # but gen_int_all and geno_cases_cc already subset them out
   # Try creating a new fun leg file that has the common variants removed
   # common_all_fun = common_all %>% filter(fun=="fun")
   # leg_fun_all_rare = subset(leg_fun, !(id %in% common_all_fun$id))
   # 
   # common_ext_fun = common_ext %>% filter(fun=="fun")
   # leg_fun_ext_rare = subset(leg_fun, !(id %in% common_ext_fun$id))
   # 
   # common_int_fun = common_int %>% filter(fun=="fun")
   # leg_fun_int_rare = subset(leg_fun, !(id %in% common_int_fun$id))
   
   # create MAC matrix for external controls
   # tbl = data.frame(a0=rowSums(geno_ext)) %>% mutate(a1=2*ncol(geno_ext)-a0)
   # tbl = data.frame(a0=geno_ext$mac) %>% mutate(a1=2*Ncc-a0) # need to use MAC
   tbl = data.frame(a0=geno_ext$mac) %>% mutate(a1=2*Ncc-a0, gene = geno_ext$gene)
   
   
   # call ProxECATv2/iECAT/SKAT once per gene
   prox2_int_genes = prox2_ext_genes = prox2_all_genes = c()
   iecat_genes = skato_int_genes = skato_ext_genes = skato_all_genes = c()
   skat_int_genes = skat_ext_genes = skat_all_genes = c()
   burden_int_genes = burden_ext_genes = burden_all_genes = c()
   
   genes = levels(droplevels(as.factor(leg$gene)))
   # g = 1
   # gene_counts = leg %>% count(gene)
   # loop through the genes
   for(g in 1:length(genes)){

     # print(paste0('current gene: ', genes[g], ' (', g, ' of ', length(genes), ')'))

     # LogProx
     # Filter data by gene
     data_int_gene = data_int %>% filter(gene==genes[g])
     data_ext_gene = data_prox %>% filter(gene==genes[g])
     data_all_gene = data_all %>% filter(gene==genes[g])
     
     # Count the number of fun and syn alleles by case status
     # need .drop param so it still creates a group even if AC is 0
     counts_data_int_gene = data_int_gene %>% count(case, fun, .drop = FALSE)
     counts_data_ext_gene = data_ext_gene %>% count(case, fun, .drop = FALSE)
     counts_data_all_gene = data_all_gene %>% count(case, fun, .drop = FALSE)
     
     # If sum of fun alleles or sum of syn alleles is < 5, mark as NA, else run LogProx
     prox2_int = ifelse(counts_data_int_gene$n[1] + counts_data_int_gene$n[3] < 5 |
                    counts_data_int_gene$n[2] + counts_data_int_gene$n[4] < 5, NA,
                  summary(glm(fun ~ case, data=data_int_gene, family="binomial"))$coefficients[2,4])

     prox2_ext = ifelse(counts_data_ext_gene$n[1] + counts_data_ext_gene$n[3] < 5 |
                        counts_data_ext_gene$n[2] + counts_data_ext_gene$n[4] < 5, NA,
                      summary(glm(fun ~ case, data=data_ext_gene, family="binomial"))$coefficients[2,4])

     prox2_all = ifelse(counts_data_all_gene$n[1] + counts_data_all_gene$n[3] < 5 |
                        counts_data_all_gene$n[2] + counts_data_all_gene$n[4] < 5, NA,
                      summary(glm(fun ~ case + group, data=data_all_gene, family="binomial"))$coefficients[2,4])
     
     # Save the LogProx p-values
     prox2_int_genes = c(prox2_int_genes, prox2_int)
     prox2_ext_genes = c(prox2_ext_genes, prox2_ext)
     prox2_all_genes = c(prox2_all_genes, prox2_all)

     ### Prepare data for iECAT and SKAT methods
     
     # Z_int = geno_int_all[which(leg_fun$gene==genes[g]), ]
     # Z_ext = geno_cases_cc[which(leg_fun$gene==genes[g]), ]

     # Z_int_all = geno_int_all[which(leg_fun_all_rare$gene==genes[g]), ]
     # Z_int = geno_cases_int[which(leg_fun_int_rare$gene==genes[g]), ]
     # Z_ext = geno_cases_cc[which(leg_fun_ext_rare$gene==genes[g]), ]
     # Z_all = geno_all[which(leg_fun_all_rare$gene==genes[g]), ]
     
     Z_int = geno_cases_int %>% filter(gene == genes[g]) %>% select(-gene)
     Z_ext = geno_cases_cc %>% filter(gene == genes[g]) %>% select(-gene)
     Z_all = geno_all %>% filter(gene == genes[g]) %>% select(-gene)
     Z_int_all = geno_int_all %>% filter(gene == genes[g]) %>% select(-gene)
     
     # subset the MAC matrix for the external controls
     # tbl_gene = tbl[which(leg_fun$gene==genes[g]), ]
     # tbl_gene = tbl[which(leg_fun_all_rare$gene==genes[g]), ]
     tbl_gene = tbl %>% filter(gene == genes[g]) %>% select(-gene)
     
     # call the iECAT-O and SKAT-O functions
     re_gene = iECAT(t(Z_int_all), obj_int, as.matrix(tbl_gene), method="optimal") # iECAT and SKAT-O internal
     skato_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="SKATO") # SKAT-O external
     skato_int_gene = SKATBinary(t(Z_int), obj_int, method="SKATO") # SKAT-O internal
     skato_all_gene = SKATBinary(t(Z_all), obj_all, method="SKATO") # SKAT-O internal+external
     
     #Add info to iecat_info
     # iecat_info <- rbind(iecat_info, c(genes[g], re_gene$param$n.marker, re_gene$param$n.marker.test, re_gene$p.value))
     
     # Save the iECAT-O and SKAT-O p-values
     iecat_genes = c(iecat_genes, re_gene$p.value)
     # iecat_skato_int_genes = c(iecat_skato_int_genes, re_gene$p.value.internal) # Check this is the same as SKAT-O, NOT THE SAME
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

   # Check counts of MAC > 2 for internal and external samples
   # mac_df = data.frame(matrix(nrow = length(genes), ncol = 3))
   # colnames(mac_df) = c("Gene", "Internal ACs > 2", "External ACs > 2")
   # for (g in 1:length(genes)) {
   #   mac_df[g, 1] = genes[g]
   #
   #   Z_int = geno_int_all[which(leg_fun$gene==genes[g]),]
   #   mac_df[g, 2] = length(which(rowSums(Z_int) > 2))
   #
   #   tbl_gene = tbl[which(leg_fun$gene==genes[g]),]
   #   mac_df[g, 3] = length(which(tbl_gene$a0 > 2))
   # }

   # for (g in 1:length(genes)) {
   #   Z_int = geno_int_all[which(leg_fun$gene==genes[g]), ]
   #   if (any(is.na(Z_int))) {
   #     print(paste("Gene ", genes[g], " contains ", sum(any(is.na(Z_int))), " NA value(s)."))
   #   }
   #   else {
   #     print(paste("Gene ", genes[g], " contains no NA values."))
   #   }
   # }

   
   # store the LogProx, iECAT-O, SKAT gene p-values
   # Each col represents a gene and each row represents a sim rep
   prox2_ext_genes_p = rbind(prox2_ext_genes_p, prox2_ext_genes)
   prox2_int_genes_p = rbind(prox2_int_genes_p, prox2_int_genes)
   prox2_all_genes_p = rbind(prox2_all_genes_p, prox2_all_genes)

   iecat_genes_p = rbind(iecat_genes_p, iecat_genes)
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

# colnames(iecat_genes_p) = c('Gene', 'n.marker', 'n.marker.test', 'p.value')
# write.table(iecat_genes_p, paste0(dir_out, "iecat_info_", folder, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)

#colnames(prox.genes.p) = colnames(prox2.genes.p) = colnames(iecat.genes.p) = colnames(skat.int.genes.p) = colnames(skat.ext.genes.p)= genes

# results = data.frame(prox.p, prox.int.p, prox.all.p, prox.all2.p, prox2.p, prox2.int.p, prox2.all.p, prox2.all2.p,
#                      iecat.p, skat.int.p, skat.ext.p, skat.all.p, skat.all2.p)

# results = data.frame(prox_ext_genes_p, prox_int_genes_p,
#                      prox2_ext_genes_p, prox2_all_genes_p, prox2_int_genes_p,
#                      iecat_genes_p, skato_int_genes_p, skato_ext_genes_p, skato_all_genes_p,
#                      skat_int_genes_p, skat_ext_genes_p, skat_all_genes_p,
#                      burden_int_genes_p, burden_ext_genes_p, burden_all_genes_p)

# Set col names to the genes
colnames(prox_ext_genes_p) = colnames(prox_int_genes_p) = genes
colnames(prox_weighted_ext_genes_p) = colnames(prox_weighted_int_genes_p) = genes
colnames(prox2_ext_genes_p) = colnames(prox2_all_genes_p) = colnames(prox2_int_genes_p) = genes
colnames(iecat_genes_p) = colnames(skato_int_genes_p) = colnames(skato_ext_genes_p) = colnames(skato_all_genes_p) = genes
colnames(skat_int_genes_p) = colnames(skat_ext_genes_p) = colnames(skat_all_genes_p) = genes
colnames(burden_int_genes_p) = colnames(burden_ext_genes_p) = colnames(burden_all_genes_p) = genes

#ProxECAT
write.table(prox_ext_genes_p, paste0(dir_out, "T1e_gene_prox_ext_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(prox_int_genes_p, paste0(dir_out, "T1e_gene_prox_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_ext_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_ext_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(prox_weighted_int_genes_p, paste0(dir_out, "T1e_gene_prox_weighted_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
# LogProx
write.table(prox2_ext_genes_p, paste0(dir_out, "T1e_gene_prox2_ext_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(prox2_all_genes_p, paste0(dir_out, "T1e_gene_prox2_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(prox2_int_genes_p, paste0(dir_out, "T1e_gene_prox2_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
# iECAT-O and SKAT-O
write.table(iecat_genes_p, paste0(dir_out, "T1e_gene_iecat_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(skato_int_genes_p, paste0(dir_out, "T1e_gene_skato_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(skato_ext_genes_p, paste0(dir_out, "T1e_gene_skato_ext_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(skato_all_genes_p, paste0(dir_out, "T1e_gene_skato_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
# SKAT
write.table(skat_int_genes_p, paste0(dir_out, "T1e_gene_skat_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(skat_ext_genes_p, paste0(dir_out, "T1e_gene_skat_ext_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(skat_all_genes_p, paste0(dir_out, "T1e_gene_skat_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
# Burden
write.table(burden_int_genes_p, paste0(dir_out, "T1e_gene_burden_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(burden_ext_genes_p, paste0(dir_out, "T1e_gene_burden_ext_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)
write.table(burden_all_genes_p, paste0(dir_out, "T1e_gene_burden_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop, "_maf", maf, ".txt"), quote=F, row.names=F, col.names=T)

#write.table(prox.p, paste0("T1e_Confounding_ProxECAT_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(prox2.p, paste0("T1e_Confounding_ProxECATv2_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(iecat.p, paste0("T1e_Confounding_iECAT_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(skat.int.p, paste0("T1e_Confounding_SKAT_internal_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(skat.ext.p, paste0("T1e_Confounding_SKAT_cc_MAF", maf, ".txt"), quote=F, row.names=F)

#write.table(prox.genes.p, paste0("T1e_Confounding_ProxECAT_genes_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(prox2.genes.p, paste0("T1e_Confounding_ProxECATv2_genes_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(iecat.genes.p, paste0("T1e_Confounding_iECAT_genes_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(skat.int.genes.p, paste0("T1e_Confounding_SKAT_internal_genes_MAF", maf, ".txt"), quote=F, row.names=F)
#write.table(skat.ext.genes.p, paste0("T1e_Confounding_SKAT_cc_genes_MAF", maf, ".txt"), quote=F, row.names=F)