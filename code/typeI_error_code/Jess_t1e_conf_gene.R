# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(proxecat)
library(SKAT)
library(iECAT)

# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")

source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")

Pop = 'NFE'
pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '100v99'
p_case = 100
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 99
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)

# dir = '/storage/math/projects/compinfo/simulations/'
# setwd(paste0(dir, 'output/20K_', Pop, '/'))

dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')

# create empty vectors to store the p-values from each replicate
# prox.p = prox.int.p = prox.all.p = c()
# prox2.p = prox2.int.p = prox2.all.p = c()
# iecat.p = skat.int.p = skat.ext.p = skat.all.p = c()
# prox_genes_p = prox2_genes_p = iecat_genes_p = skat_int_genes_p = skat_ext_genes_p = c()

prox_ext_genes_p = prox_int_genes_p = c() #proxECAT
prox2_ext_genes_p = prox2_int_genes_p = prox2_all_genes_p = c() #LogProx
iecat_genes_p = skato_ext_genes_p = skato_int_genes_p = skato_all_genes_p = c() #iECAT-O and SKAT-O
skat_ext_genes_p = skat_int_genes_p = skat_all_genes_p = c() #SKAT
burden_ext_genes_p = burden_int_genes_p = burden_all_genes_p = c() #Burden


# loop through the simulation replicates
set.seed(1) 
i=1
for (i in 1:100){
  
   # read in the legend file
   # leg = read_leg_homo(dir_leg, Pop, i)
   leg = read.table(paste0(dir_leg, 'chr19.block37.', Pop, '.sim', i, '.', p_case, 'fun.', p_case, 'syn.legend'), header=T, sep='\t') #RAREsim v2.1.1 pruning only
   leg$row = 1:nrow(leg)
   
   # read in the haplotype files
   hap_cases = read_hap_homo(dir_in, Pop, i, "cases", p_case_fun, p_case_syn)
   hap_int = read_hap_homo(dir_in, Pop, i, "internal.controls", p_int_fun, p_int_syn)
   hap_cc = read_hap_homo(dir_in, Pop, i, "common.controls", p_cc_fun, p_cc_syn)

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
   counts_ext_gene = data_prox %>% count(gene, case, fun)
   counts_ext_wide = tidyr::pivot_wider(counts_ext_gene, names_from=c(case, fun), values_from=n,
                                        values_fill=0, names_sep="_")
   
   counts_int_gene = data_int %>% count(gene, case, fun)
   counts_int_wide = tidyr::pivot_wider(counts_int_gene, names_from=c(case, fun), values_from=n,
                                        values_fill=0, names_sep="_")
   # counts_all = colSums(counts_wide[,-1])
   # names(counts_all) = colnames(counts_wide)[-1]

   # counts_all = data_all %>% count(case, fun)
   # counts_prox = data_prox %>% count(case, fun)
   # counts_int = data_int %>% count(case, fun)

   # call proxECAT
   #prox = proxecat(counts.all['case.fun'], counts.all['case.syn'], counts.all['control.fun'], counts.all['control.syn'])
   # prox.genes = proxecat(counts_wide$case_fun, counts_wide$case_syn, counts_wide$control_fun, counts_wide$control_syn)

   # prox = proxecat(counts.prox$n[1], counts.prox$n[2], counts.prox$n[3], counts.prox$n[4])
   # prox.all = proxecat(counts.all$n[1], counts.all$n[2], counts.all$n[3], counts.all$n[4])
   # prox.int = proxecat(counts.int$n[1], counts.int$n[2], counts.int$n[3], counts.int$n[4])

   # store proxECAT p-values
   # prox.p = rbind(prox.p, prox$p.value)
   # prox.all.p = rbind(prox.all.p, prox.all$p.value)
   # prox.int.p = rbind(prox.int.p, prox.int$p.value)
   # prox.genes.p = rbind(prox.genes.p, prox.genes$p.value)
   
   # create case/control phenotype matrices for iECAT/SKAT
   pheno_int = rep(0, (ncol(geno_cases) + ncol(geno_int)))
   pheno_int[1:ncol(geno_cases)] = 1

   pheno_ext = rep(0, (ncol(geno_cases) + ncol(geno_cc)))
   pheno_ext[1:ncol(geno_cases)] = 1

   pheno_all = rep(0, (ncol(geno_cases) + ncol(geno_int) + ncol(geno_cc)))
   pheno_all[1:ncol(geno_cases)] = 1
   
   # subset the synonymous variants from the legend file
   leg_syn = leg %>% filter(fun=="syn")
   leg_fun = leg %>% filter(fun=="fun")
   
   # create combined genotype matrices
   geno_int_all = cbind(geno_cases, geno_int)[-union(leg_syn$row, common_all$row),] # for iECAT
   geno_cases_cc = cbind(geno_cases, geno_cc)[-union(leg_syn$row, common_ext$row),] # SKAT (just cases & external controls)
   # geno_ext = geno_cc[-union(leg_syn$row, common_all$row),] #will need to change geno_cc to counts_cc for admixed
   geno_ext = count_cc[-union(leg_syn$row, common_all$row),] # need MAC for tbl object for iECAT
   geno_all = cbind(geno_cases, geno_int, geno_cc)[-union(leg_syn$row, common_all$row),] # SKAT (all)

   # null model object
   obj_int = SKAT_Null_Model(as.numeric(pheno_int) ~ 1, out_type="D") # D-dichotomous
   obj_ext = SKAT_Null_Model(as.numeric(pheno_ext) ~ 1, out_type="D") # D-dichotomous
   obj_all = SKAT_Null_Model(as.numeric(pheno_all) ~ 1, out_type="D") # D-dichotomous

   # fit the ProxECATv2 model
   #glm.prox = glm(fun ~ case + group, data=data.all, family="binomial") 
   # glm.prox = glm(fun ~ case, data=data.prox, family="binomial") 
   # glm.all.prox = glm(fun ~ case + group, data=data.all, family="binomial") 
   # glm.int.prox = glm(fun ~ case, data=data.int, family="binomial") 
   
   # save the p-value for case/control status
   # p.prox = summary(glm.prox)$coefficients[2,4]
   # p.prox.all = summary(glm.all.prox)$coefficients[2,4]
   # p.prox.int = summary(glm.int.prox)$coefficients[2,4]
   
   # create MAC matrix for external controls
   # tbl = data.frame(a0=rowSums(geno_ext)) %>% mutate(a1=2*ncol(geno_ext)-a0)
   tbl = data.frame(a0=geno_ext$mac) %>% mutate(a1=2*Ncc-a0) # need to use MAC
   
   # call the iECAT function
   # re = iECAT(t(geno.int.all), obj.int, as.matrix(tbl), method="optimal")
   # re.skat = SKATBinary(t(geno.cases.cc), obj.ext, method="SKATO") # SKAT-O based on the unified approach
   # re.all = SKATBinary(t(geno.all), obj.all, method="SKATO") # SKAT-O based on the unified approach
   
   # save the p-values
   # prox2.p = c(prox2.p, p.prox)
   # prox2.all.p = c(prox2.all.p, p.prox.all)
   # prox2.int.p = c(prox2.int.p, p.prox.int)
   # iecat.p = c(iecat.p, re$p.value)
   # skat.int.p = c(skat.int.p, re$p.value.internal)
   # skat.ext.p = c(skat.ext.p, re.skat$p.value)
   # skat.all.p = c(skat.all.p, re.all$p.value)
   
   # call ProxECATv2/iECAT/SKAT once per gene
   prox_ext_genes = prox_int_genes = c()
   prox2_int_genes = prox2_ext_genes = prox2_all_genes = c()
   iecat_genes = skato_int_genes = skato_ext_genes = skato_all_genes = c()
   skat_int_genes = skat_ext_genes = skat_all_genes = c()
   burden_int_genes = burden_ext_genes = burden_all_genes =c()
   
   genes = levels(droplevels(as.factor(leg$gene)))
   g = 1
   gene_counts = leg %>% count(gene)
   # loop through the genes
   for(g in 1:length(genes)){ 

     print(paste0('current gene: ', genes[g]))
     
     # subset the data by gene
     # ProxECAT
     counts_ext_wide_gene = counts_ext_wide %>% filter(gene==genes[g])
     counts_int_wide_gene = counts_int_wide %>% filter(gene==genes[g])
     
     # LogProx
     data_all_gene = data_all %>% filter(gene==genes[g])
     data_ext_gene = data_prox %>% filter(gene==genes[g])
     data_int_gene = data_int %>% filter(gene==genes[g])

     # return NA if there are no minor alleles in any of the groups
     if(summary(data_all_gene$group)[[2]]==0 | summary(data_all_gene$group)[[1]]==0 | #internal or external
        summary(data_all_gene$case)[[2]]==0 | summary(data_all_gene$case)[[1]]==0 |   #control or case
        summary(data_all_gene$fun)[[2]]==0 | summary(data_all_gene$fun)[[1]]==0){     #syn or fun

       prox2_int_genes = c(prox2_int_genes, NA)
       prox2_ext_genes = c(prox2_ext_genes, NA)
       prox2_all_genes = c(prox2_all_genes, NA)
       
       iecat_genes = c(iecat_genes, NA)
       skato_int_genes = c(skato_int_genes, NA)
       skato_ext_genes = c(skato_ext_genes, NA)
       skato_all_genes = c(skato_all_genes, NA)
       
       skat_int_genes = c(skat_int_genes, NA)
       skat_ext_genes = c(skat_ext_genes, NA)
       skat_all_genes = c(skat_all_genes, NA)
       burden_int_genes = c(burden_int_genes, NA)
       burden_ext_genes = c(burden_ext_genes, NA)
       burden_all_genes = c(burden_all_genes, NA)

     } else {

        # Run proxECAT
        prox_ext_gene = proxecat(counts_ext_wide_gene$case_fun, counts_ext_wide_gene$case_syn,
                                 counts_ext_wide_gene$control_fun, counts_ext_wide_gene$control_syn)
        prox_int_gene = proxecat(counts_int_wide_gene$case_fun, counts_int_wide_gene$case_syn, 
                                 counts_int_wide_gene$control_fun, counts_int_wide_gene$control_syn)
        
        # Save the proxECAT p-values
        prox_ext_genes = c(prox_ext_genes, prox_ext_gene$p.value)
        prox_int_genes = c(prox_int_genes, prox_int_gene$p.value)
       
        # fit the ProxECATv2 model
        glm_ext_prox = glm(fun ~ case + group, data=data_ext_gene, family="binomial")
        glm_int_prox = glm(fun ~ case + group, data=data_int_gene, family="binomial")
        glm_all_prox = glm(fun ~ case + group, data=data_all_gene, family="binomial")

        # Save the LogProx p-values
        prox2_ext_genes = c(prox2_ext_genes, summary(glm_ext_prox)$coefficients[2,4])
        prox2_int_genes = c(prox2_int_genes, summary(glm_int_prox)$coefficients[2,4])
        prox2_all_genes = c(prox2_all_genes, summary(glm_all_prox)$coefficients[2,4])

        # subset the genotype matrices
        # LAST FEW ROWS ARE NA FOR THE FIRST GENE, NEED TO FIGURE OUT WHY
        # I think it's because leg_fun still has the common variants in it
        # but gen_int_all and geno_cases_cc already subset them out
        # Try creating a new fun leg file that has the common variants removed
        common_all_fun = common_all %>% filter(fun=="fun")
        leg_fun_all_rare = subset(leg_fun, !(id %in% common_all_fun$id))
        
        common_ext_fun = common_ext %>% filter(fun=="fun")
        leg_fun_ext_rare = subset(leg_fun, !(id %in% common_ext_fun$id))
        # Z_int = geno_int_all[which(leg_fun$gene==genes[g]), ]
        # Z_ext = geno_cases_cc[which(leg_fun$gene==genes[g]), ]
        Z_int = geno_int_all[which(leg_fun_all_rare$gene==genes[g]), ]
        Z_ext = geno_cases_cc[which(leg_fun_ext_rare$gene==genes[g]), ]
        Z_all = geno_all[which(leg_fun_all_rare$gene==genes[g]), ]
        
        print(paste0("Length of Z_int: ", nrow(Z_int)))
        print(paste0("Length of Z_ext: ", nrow(Z_int)))

        # subset the MAC matrix for the external controls
        # tbl_gene = tbl[which(leg_fun$gene==genes[g]), ]
        tbl_gene = tbl[which(leg_fun_all_rare$gene==genes[g]), ]

        # call the iECAT-O and SKAT-O functions
        re_gene = iECAT(t(Z_int), obj_int, as.matrix(tbl_gene), method="optimal") # iECAT and SKAT-O internal
        skato_ext_gene = SKATBinary(t(Z_ext), obj_ext, method="SKATO") # SKAT-O external
        skato_all_gene = SKATBinary(t(Z_all), obj_all, method="SKATO") # SKAT-O internal+external
        
        # Save the iECAT-O and SKAT-O p-values
        iecat_genes = c(iecat_genes, re_gene$p.value)
        skato_int_genes = c(skato_int_genes, re_gene$p.value.internal)
        skato_ext_genes = c(skato_ext_genes, skato_ext_gene$p.value)
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

   # store the ProxECAT, LogProx, iECAT-O, SKAT gene p-values
   # Each col represents a gene and each row represents a sim rep
   prox_ext_genes_p = rbind(prox_ext_genes_p, prox_ext_genes)
   prox_int_genes_p = rbind(prox_int_genes_p, prox_int_genes)
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

#colnames(prox.genes.p) = colnames(prox2.genes.p) = colnames(iecat.genes.p) = colnames(skat.int.genes.p) = colnames(skat.ext.genes.p)= genes

results = data.frame(prox.p, prox.int.p, prox.all.p, prox.all2.p, prox2.p, prox2.int.p, prox2.all.p, prox2.all2.p,
                     iecat.p, skat.int.p, skat.ext.p, skat.all.p, skat.all2.p)

setwd(paste0(dir, 'results/20K_', Pop, '/'))

write.table(results, paste0("T1e_Confounding_", Pop, "_maf", maf, "_", end, ".txt"), quote=F, row.names=F)

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