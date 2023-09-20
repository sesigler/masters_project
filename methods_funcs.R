############################################################################## 
# This file contains the functions necessary to format the data, run the statistical 
# test, and extract the p-value for each of the 3 methods (proxECAT, LogProx, iECAT-O)
#
# Variable legend:
# leg: the legend file
# counts.<DATASET>: dataframe containing the allele counts, MACs, and MAFs of
#                   DATASET'S genoytpe file     
# common: the common variants that are common in at least one of the relevant datasets
# common.<ext/all>: common variants for either cases+cc or cases+int+cc
# adj: Boolean value indicating whether or not the adjusted values are being used
# geno.<DATASET>: genotype matrix for DATASET
# Ncc: Number of individuals in the common controls
##############################################################################

### Determine the number of rare fun and syn minor alleles in a dataset
rare_counts = function(counts, leg.fun, leg.syn, maf){
  
  fun.counts = counts[leg.fun$row, ]
  rare.fun = which(fun.counts$maf <= maf)
  out.fun = sum(fun.counts[rare.fun, ]$mac)
  
  syn.counts = counts[leg.syn$row, ]
  rare.syn = which(syn.counts$maf <= maf)
  out.syn = sum(syn.counts[rare.syn, ]$mac)
  
  out = c(out.fun, out.syn)
  
  return(out)
}


# Function for formatting data and running statistical test for ProxECAT
prox_data_prep = function(leg.fun, leg.syn, counts.cases, counts.ctrl, maf) {
  
  counts.prox = c()
  
  case.fun = rare_counts(counts.cases, leg.fun, leg.syn, maf)
  counts.prox = c(counts.prox, c(case.fun[1], case.fun[2]))
  
  ctrl.fun = rare_counts(counts.ctrl, leg.fun, leg.syn, maf)
  counts.prox = c(counts.prox, c(ctrl.fun[1], ctrl.fun[2]))
  
  # Run proxECAT
  prox = proxecat(counts.prox[1], counts.prox[2], counts.prox[3], counts.prox[4])
  
  # return p-value
  return(prox$p.value)
}

# # Function for formatting data and running statistical test for ProxECAT
# prox_data_prep = function(leg, counts.cases, counts.cc, common, adj) {
#   
#   # convert genotypes into long format for ProxECAT v2
#   data.cases = make_long(counts.cases, leg, "case", "int")
#   if (adj){
#     data.cc = make_long_adj(counts.cc, leg, "control", "ext") #doesn't have count column
#   }
#   else{
#     data.cc = make_long(counts.cc, leg, "control", "ext")
#   }
#   
#   # combine the data together AND REMOVE COMMON VARIANTS
#   data.prox = data.frame(lapply(rbind(data.cases, data.cc), factor)) %>% 
#     filter(!(id %in% common$id))
#   
#   # getting overall counts for functional & case status
#   # data for proxECAT method
#   counts.prox = data.prox %>% count(case, fun)
#   
#   # Run proxECAT
#   prox = proxecat(counts.prox$n[1], counts.prox$n[2], counts.prox$n[3], counts.prox$n[4])
#   
#   # return p-value
#   return(prox$p)
# }

# Function for formatting data and running statistical test for ProxECAT
# for testing cases vs internal controls only
# prox_int_prep = function(leg, counts.cases, counts.int, common) {
#   
#   # convert genotypes into long format for ProxECAT v2
#   data.cases = make_long(counts.cases, leg, "case", "int")
#   data.int = make_long(counts.int, leg, "control", "int")
#  
#   # combine the data together AND REMOVE COMMON VARIANTS
#   data.prox = data.frame(lapply(rbind(data.cases, data.int), factor)) %>% 
#     filter(!(id %in% common$id))
#   
#   # getting overall counts for functional & case status
#   # data for proxECAT method
#   counts.prox = data.prox %>% count(case, fun)
#   
#   # Run proxECAT
#   prox = proxecat(counts.prox$n[1], counts.prox$n[2], counts.prox$n[3], counts.prox$n[4])
#   
#   # return p-value
#   return(prox$p)
# }


# Function for formatting data and running statistical test for LogProx
logprox_data_prep = function(leg, counts.cases, counts.int, counts.cc, common.ext, common.all, adj) {
  
  # convert genotypes into long format for ProxECAT v2
  data.cases = make_long(counts.cases, leg, "case", "int")
  data.int = make_long(counts.int, leg, "control", "int")
  
  if (adj) {
    data.cc = make_long_adj(counts.cc, leg, "control", "ext") #doesn't have count column
  }
  else {
    data.cc = make_long(counts.cc, leg, "control", "ext")
  }
  
  # combine the data together AND REMOVE COMMON VARIANTS
  data.prox = data.frame(lapply(rbind(data.cases, data.cc), factor)) %>% 
    filter(!(id %in% common.ext$id))
  
  data.all = data.frame(lapply(rbind(data.cases, data.int, data.cc), factor)) %>% 
    filter(!(id %in% common.all$id))
  
  # fit the ProxECATv2 model
  glm.prox = glm(fun ~ case, data=data.prox, family="binomial") 
  glm.all.prox = glm(fun ~ case + group, data=data.all, family="binomial")
  
  # save the p-value for case/control status
  p.prox = summary(glm.prox)$coefficients[2,4]
  p.prox.all = summary(glm.all.prox)$coefficients[2,4]
  
  out <- list(p.prox, p.prox.all)
  
  return(out) 
  
}

# Function for formatting data and running statistical test for LogProx
# for testing cases vs internal controls only
logprox_int_prep = function(leg, counts.cases, counts.int, common) {
  
  # convert genotypes into long format for ProxECAT v2
  data.cases = make_long(counts.cases, leg, "case", "int")
  data.int = make_long(counts.int, leg, "control", "int")
  
  # combine the data together AND REMOVE COMMON VARIANTS
  data.int = data.frame(lapply(rbind(data.cases, data.int), factor)) %>% 
    filter(!(id %in% common$id))
  
  # fit the ProxECATv2 model
  glm.int = glm(fun ~ case, data=data.int, family="binomial") 
  
  # save the p-value for case/control status
  p.int = summary(glm.int)$coefficients[2,4]
  
  return(p.int) 
  
}


# Function for formatting data and running statistical test for iECAT-O
iecat_data_prep = function(geno.cases, geno.int, leg, common, counts.cc, Ncc) {
  
  # create case/control phenotype matrices for iECAT/SKAT
  pheno.int = rep(0, (ncol(geno.cases) + ncol(geno.int))) 
  pheno.int[1:ncol(geno.cases)] = 1
  
  # subset the synonymous variants from the legend file
  leg.syn = leg %>% filter(fun=="syn")
  
  # create combined genotype matrices and remove the synonymous & common variants
  geno.int.all = cbind(geno.cases, geno.int)[-union(leg.syn$row, common$row),]
  
  geno.ext = counts.cc[-union(leg.syn$row, common$row),]
  
  # null model object
  # distinguishes between cases and internal controls
  obj.int = SKAT_Null_Model(as.numeric(pheno.int) ~ 1, out_type="D") # D-dichotomous
  
  # create MAC matrix for external controls
  # geno.ext.adj needs to be adj MAC counts but only for FUNC variants and RARE variants, so subset those
  # a0 is MAC. a1 is Major ACs, so it's just ncol(hap)-a0
  # a0+a1 should sum to 20000 bc full hap file
  tbl = data.frame(a0=geno.ext$mac) %>% mutate(a1=2*Ncc-a0)
  
  # call the iECAT function
  re = iECAT(t(geno.int.all), obj.int, as.matrix(tbl), method="optimal")
  
  # extract the p-values (iECAT-O and SKAT-O internal)
  out = list(re$p.value, re$p.value.internal)
  
  return(out)
}


# Function for formatting data and running statistical test for SKAT-O
skat_data_prep = function(geno.cases, geno.int, geno.cc, leg, common.ext, common.all) {
  
  # create case/control phenotype matrices for SKAT
  pheno.ext = rep(0, (ncol(geno.cases) + ncol(geno.cc))) 
  pheno.ext[1:ncol(geno.cases)] = 1
  
  pheno.all = rep(0, (ncol(geno.cases) + ncol(geno.int) + ncol(geno.cc))) 
  pheno.all[1:ncol(geno.cases)] = 1
  
  # subset the synonymous variants from the legend file
  leg.syn = leg %>% filter(fun=="syn")
  
  # create combined genotype matrices and remove the synonymous & common variants
  geno.cases.cc = cbind(geno.cases, geno.cc)[-union(leg.syn$row, common.ext$row),] # SKAT (just cases & external controls)
  geno.all = cbind(geno.cases, geno.int, geno.cc)[-union(leg.syn$row, common.all$row),] # SKAT (all)
  
  # null model object
  obj.ext = SKAT_Null_Model(as.numeric(pheno.ext) ~ 1, out_type="D") # D-dichotomous
  obj.all = SKAT_Null_Model(as.numeric(pheno.all) ~ 1, out_type="D") # D-dichotomous
  
  # call the iECAT function
  re.skat = SKATBinary(t(geno.cases.cc), obj.ext, method="SKATO") # SKAT-O based on the unified approach
  re.all = SKATBinary(t(geno.all), obj.all, method="SKATO") # SKAT-O based on the unified approach   
  
  # extract the p-valus
  out = list(re.skat$p.value, re.all$p.value)
  
  return(out)
}