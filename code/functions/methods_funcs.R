############################################################################## 
# This file contains the functions necessary to format the data, run the statistical 
# test, and extract the p-value for each of the 3 methods (proxECAT, LogProx, iECAT-O)
##############################################################################

#' prox_gene_data_prep
#' 
#' @description
#' Function for formatting data and running statistical test for ProxECAT by gene
#' 
#' @param count.cases the dataframe of ACs and AFs for each variant in the cases
#' @param count.controls the dataframe of ACs and AFs for each variant in the controls
#' @param leg the legend file
#' @param common the legend file filtered down to variants that are common in 
#'               either the cases or controls
#' 
#' @return a dataframe containing the ACs by gene, functional, and case status 
#'         as well as the p-values for proxECAT and proxECAT-weighted

prox_gene_data_prep = function(count.cases, count.controls, leg, common) {
  
  # data.cases = count.cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case", group="int") %>% 
  #   filter(mac!=0) %>% select(-count)
  data.cases = count.cases %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="case") %>% 
    filter(ac!=0)
  
  data.controls = count.controls %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case="control") %>% 
    filter(ac!=0)

  data.prox = data.frame(rbind(data.cases, data.controls)) %>% mutate(across(all_of(c("id", "gene", "fun", "case")), as.factor)) %>% 
    filter(!(id %in% common$id))
  
  counts_gene = data.prox %>% group_by(gene, case, fun) %>% summarise(n = sum(ac))

  counts_wide = tidyr::pivot_wider(counts_gene, names_from=c(case, fun), values_from=n,
                                   values_fill=0, names_sep="_")
  
  # Calculate the ratios
  counts_wide2 = counts_wide %>% mutate(case_ratio = case_fun/case_syn,
                                       control_ratio = control_fun/control_syn)

  # Calculate medians
  # Set na.rm to TRUE to avoid median ratio being NA if there's a divide by 0 problem
  median_case_ratio = median(counts_wide2$case_ratio, na.rm = TRUE)
  median_control_ratio = median(counts_wide2$control_ratio, na.rm = TRUE)

  # Calculate the weighted values and the p-values
  counts_wide3 = counts_wide2 %>% mutate(case_fun_w = case_fun / median_case_ratio,
                                         control_fun_w = control_fun / median_control_ratio) %>%
    mutate(prox = ifelse(case_fun + control_fun < 5 | case_syn + control_syn < 5, NA,
                         proxecat(case_fun, case_syn, control_fun, control_syn)$p.value),
           prox_w = ifelse(case_fun_w + control_fun_w < 5 | case_syn + control_syn < 5, NA,
                           proxecat(case_fun_w, case_syn, control_fun_w, control_syn)$p.value))
  
  return(counts_wide3)
}

#' format_logprox_data
#' 
#' @description
#' Function to format the data used for LogProx
#' 
#' @param leg the legend file
#' @param count.case the dataframe of ACs and AFs for each variant in the cases
#' @param count.control the dataframe of ACs and AFs for each variant in the controls
#' @param control_type a string denoting whether count.control is int (internal) 
#'                     or ext (external)
#' @param count.control2 the dataframe of ACs and AFs for each variant in the 
#'                       common controls if internal and external data are used, 
#'                       default value is NULL
#' @param common the legend file filtered down to variants that are common in 
#'               either the cases or controls
#' @param data.all boolean value denoting if both internal and external controls 
#'                 are being used
#' 
#' @return a dataframe with id, gene, fun, case, and group as columns and a row 
#'         repeated AC times for each variant in the cases and controls

format_logprox_data = function(leg, count.case, count.control, control_type, count.control2=NULL, common, data.all=FALSE) {
  
  # Convert datasets to long format
  data.case = make_long(count.case, leg, "case", "int")
  data.control = make_long(count.control, leg, "control", control_type)
  
  if (data.all) {
    # Make long common control data
    data.control2 = make_long(count.control2, leg, "control", "ext")
    
    # cbind all data and filter out common variants
    data.out = data.frame(lapply(rbind(data.case, data.control, data.control2), factor)) %>%
      filter(!(id %in% common$id))
    
    return(data.out)
  }
  
  # cbind data and filter out common variants
  data.out = data.frame(lapply(rbind(data.case, data.control), factor)) %>%
    filter(!(id %in% common$id))
  
  return(data.out)
}

#' logprox_gene_data_prep
#' 
#' @description
#' Function to run the statistical test for LogProx by gene
#' 
#' @param data.prox2 the formatted datafrom produced by format_logprox_data
#' @param current.gene a string of the gene to run LogProx on
#' @param data.all boolean value denoting if both internal and external controls are being used
#' 
#' @return the p-value returned by LogProx

# Function for formatting data and running statistical test for LogProx
logprox_gene_data_prep = function(data.prox2, current.gene, data.all) {
  
  # LogProx
  # Filter data by gene
  data.gene = data.prox2 %>% filter(gene==current.gene)
  
  # Count the number of fun and syn alleles by case status
  # need .drop param so it still creates a group even if AC is 0
  counts.data.gene = data.gene %>% count(case, fun, .drop = FALSE)
  
  # If sum of fun alleles or sum of syn alleles is < 5, mark as NA, else run LogProx
  if (data.all) {
    prox2 = ifelse(counts.data.gene$n[1] + counts.data.gene$n[3] < 5 | 
                     counts.data.gene$n[2] + counts.data.gene$n[4] < 5, NA, 
                   summary(glm(fun ~ case + group, data=data.gene, family="binomial"))$coefficients[2,4])
    return(prox2)
  } 
  
  prox2 = ifelse(counts.data.gene$n[1] + counts.data.gene$n[3] < 5 |
                   counts.data.gene$n[2] + counts.data.gene$n[4] < 5, NA,
                 summary(glm(fun ~ case, data=data.gene, family="binomial"))$coefficients[2,4])

  
  return(prox2) 
  
}


# Function for formatting data and running statistical test for LogProx
# logprox_data_prep = function(leg, counts.cases, counts.int, counts.cc, common.ext, common.all, adj) {
#   
#   # convert genotypes into long format for ProxECAT v2
#   data.cases = make_long(counts.cases, leg, "case", "int")
#   data.int = make_long(counts.int, leg, "control", "int")
#   
#   if (adj) {
#     data.cc = make_long_adj(counts.cc, leg, "control", "ext") #doesn't have count column
#   }
#   else {
#     data.cc = make_long(counts.cc, leg, "control", "ext")
#   }
#   
#   # combine the data together AND REMOVE COMMON VARIANTS
#   data.prox = data.frame(lapply(rbind(data.cases, data.cc), factor)) %>% 
#     filter(!(id %in% common.ext$id))
#   
#   data.all = data.frame(lapply(rbind(data.cases, data.int, data.cc), factor)) %>% 
#     filter(!(id %in% common.all$id))
#   
#   # fit the ProxECATv2 model
#   glm.prox = glm(fun ~ case, data=data.prox, family="binomial") 
#   glm.all.prox = glm(fun ~ case + group, data=data.all, family="binomial")
#   
#   # save the p-value for case/control status
#   p.prox = summary(glm.prox)$coefficients[2,4]
#   p.prox.all = summary(glm.all.prox)$coefficients[2,4]
#   
#   out <- list(p.prox, p.prox.all)
#   
#   return(out) 
#   
# }

# Function for formatting data and running statistical test for LogProx
# for testing cases vs internal controls only
# logprox_int_prep = function(leg, counts.cases, counts.int, common) {
#   
#   # convert genotypes into long format for ProxECAT v2
#   data.cases = make_long(counts.cases, leg, "case", "int")
#   data.int = make_long(counts.int, leg, "control", "int")
#   
#   # combine the data together AND REMOVE COMMON VARIANTS
#   data.int = data.frame(lapply(rbind(data.cases, data.int), factor)) %>% 
#     filter(!(id %in% common$id))
#   
#   # fit the ProxECATv2 model
#   glm.int = glm(fun ~ case, data=data.int, family="binomial") 
#   
#   # save the p-value for case/control status
#   p.int = summary(glm.int)$coefficients[2,4]
#   
#   return(p.int) 
#   
# }

# Function for formatting data and running statistical test for iECAT-O
# Changed so that I can filter by either fun or syn variants
# NOTE: Whichever leg I input is what gets filtered out
# iecat_data_prep = function(geno.cases, geno.int, leg, common, counts.cc, Ncc) {
#   
#   # create case/control phenotype matrices for iECAT/SKAT
#   pheno.int = rep(0, (ncol(geno.cases) + ncol(geno.int))) 
#   pheno.int[1:ncol(geno.cases)] = 1
#   
#   # subset the synonymous variants from the legend file
#   # leg.syn = leg %>% filter(fun=="syn")
#   
#   # create combined genotype matrices and remove the synonymous & common variants
#   # geno.int.all = cbind(geno.cases, geno.int)[-union(leg.syn$row, common$row),]
#   # geno.ext = counts.cc[-union(leg.syn$row, common$row),]
#   geno.int.all = cbind(geno.cases, geno.int)[-union(leg$row, common$row),]
#   geno.ext = counts.cc[-union(leg$row, common$row),]
#   
#   # null model object
#   # distinguishes between cases and internal controls
#   obj.int = SKAT_Null_Model(as.numeric(pheno.int) ~ 1, out_type="D") # D-dichotomous
#   
#   # create MAC matrix for external controls
#   # geno.ext.adj needs to be adj MAC counts but only for FUNC variants and RARE variants, so subset those
#   # a0 is MAC. a1 is Major ACs, so it's just ncol(hap)-a0
#   # a0+a1 should sum to 20000 bc full hap file
#   tbl = data.frame(a0=geno.ext$mac) %>% mutate(a1=2*Ncc-a0)
#   
#   # call the iECAT function
#   re = iECAT(t(geno.int.all), obj.int, as.matrix(tbl), method="optimal")
#   
#   # extract the p-values (iECAT-O and SKAT-O internal)
#   out = list(re$p.value, re$p.value.internal)
#   
#   return(out)
# }


# Function for formatting data and running statistical test for SKAT-O
# Changed so that I can filter by either fun or syn variants
# NOTE: Whichever leg I input is what gets filtered out
# skato_data_prep = function(geno.cases, geno.int, geno.cc, leg, common.ext, common.all) {
#   
#   # create case/control phenotype matrices for SKAT
#   pheno.ext = rep(0, (ncol(geno.cases) + ncol(geno.cc))) 
#   pheno.ext[1:ncol(geno.cases)] = 1
#   
#   pheno.all = rep(0, (ncol(geno.cases) + ncol(geno.int) + ncol(geno.cc))) 
#   pheno.all[1:ncol(geno.cases)] = 1
#   
#   # subset the synonymous variants from the legend file
#   # leg.syn = leg %>% filter(fun=="syn")
#   
#   # create combined genotype matrices and remove the synonymous & common variants
#   # geno.cases.cc = cbind(geno.cases, geno.cc)[-union(leg.syn$row, common.ext$row),] # SKAT (just cases & external controls)
#   # geno.all = cbind(geno.cases, geno.int, geno.cc)[-union(leg.syn$row, common.all$row),] # SKAT (all)
#   geno.cases.cc = cbind(geno.cases, geno.cc)[-union(leg$row, common.ext$row),] # SKAT (just cases & external controls)
#   geno.all = cbind(geno.cases, geno.int, geno.cc)[-union(leg$row, common.all$row),] # SKAT (all)
#   
#   # null model object
#   obj.ext = SKAT_Null_Model(as.numeric(pheno.ext) ~ 1, out_type="D") # D-dichotomous
#   obj.all = SKAT_Null_Model(as.numeric(pheno.all) ~ 1, out_type="D") # D-dichotomous
#   
#   # call the SKAT function
#   re.skat = SKATBinary(t(geno.cases.cc), obj.ext, method="SKATO") # SKAT-O based on the unified approach
#   re.all = SKATBinary(t(geno.all), obj.all, method="SKATO") # SKAT-O based on the unified approach   
#   
#   # extract the p-valus
#   out = list(re.skat$p.value, re.all$p.value)
#   
#   return(out)
# }

# Function for formatting data and running statistical test for SKAT or Burden tests
# Changed so that I can filter by either fun or syn variants
# NOTE: Whichever leg I input is what gets filtered out
# skat_data_prep = function(geno.cases, geno.int, geno.cc, leg, common.int, common.ext, common.all, meth) {
#   
#   # create case/control phenotype matrices for SKAT
#   pheno.int = rep(0, (ncol(geno.cases) + ncol(geno.int))) 
#   pheno.int[1:ncol(geno.cases)] = 1 
#   
#   pheno.ext = rep(0, (ncol(geno.cases) + ncol(geno.cc))) 
#   pheno.ext[1:ncol(geno.cases)] = 1
#   
#   pheno.all = rep(0, (ncol(geno.cases) + ncol(geno.int) + ncol(geno.cc))) 
#   pheno.all[1:ncol(geno.cases)] = 1
#   
#   # subset the synonymous variants from the legend file
#   # leg.syn = leg %>% filter(fun=="syn")
#   
#   # create combined genotype matrices and remove the synonymous & common variants
#   # geno.cases.cc = cbind(geno.cases, geno.cc)[-union(leg.syn$row, common.ext$row),] # SKAT (just cases & external controls)
#   # geno.all = cbind(geno.cases, geno.int, geno.cc)[-union(leg.syn$row, common.all$row),] # SKAT (all)
#   geno.cases.int = cbind(geno.cases, geno.int)[-union(leg$row, common.int$row),] # SKAT internal 
#   geno.cases.cc = cbind(geno.cases, geno.cc)[-union(leg$row, common.ext$row),] # SKAT external
#   geno.all = cbind(geno.cases, geno.int, geno.cc)[-union(leg$row, common.all$row),] # SKAT (all)
#   
#   # null model object
#   obj.int = SKAT_Null_Model(as.numeric(pheno.int) ~ 1, out_type="D") # D-dichotomous
#   obj.ext = SKAT_Null_Model(as.numeric(pheno.ext) ~ 1, out_type="D") # D-dichotomous
#   obj.all = SKAT_Null_Model(as.numeric(pheno.all) ~ 1, out_type="D") # D-dichotomous
#   
#   # call the SKAT function
#   re.int = SKATBinary(t(geno.cases.int), obj.int, method=meth) # internal
#   re.ext = SKATBinary(t(geno.cases.cc), obj.ext, method=meth) # external
#   re.all = SKATBinary(t(geno.all), obj.all, method=meth) # all
#   
#   
#   # extract the p-valus
#   out = list(re.int$p.value, re.ext$p.value, re.all$p.value)
#   
#   return(out)
# }