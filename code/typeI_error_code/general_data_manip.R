############################################################################## 
# This file contains the functions necessary to do general data manipulation
# on the files necessary for performing type I error and power calculations for 
# several different rare variant association tests
##############################################################################

#' make_geno
#' 
#' @description
#' function to convert the haplotypes into a genotype matrix
#' 
#' @param hap A matrix of 0s and 1s where the rows represent a variant and the columns represent a single haplotype. Every 2 columns correspond to one individual
#' 
#' @return A genotype matrix where each row corresponds to an individual with values being 0, 1, or 2 reference alleles

make_geno = function(hap) {
  
  # create an empty genotype matrix
  geno = matrix(0, nrow(hap), ncol(hap)/2)
  
  # sum up the number of alleles in adjacent haplotypes (2 haplotypes per person)
  for (j in 1:(ncol(hap)/2)) {
    geno[,j] = hap[,2*j] + hap[,2*j-1]
  }
  geno = as.data.frame(geno)
  
  return(geno)
}

#' make_long
#' 
#' @description
#' function to create a dataframe with a line for each variant instead of just counts (necessary for LogProx)
#' 
#' @param counts dataframe containing the ACs and AFs of a genotype file
#' @param leg the legend file
#' @param case a string indicating if the counts are "cases" or "controls"
#' @param group a string indicating if the counts are from an "int" (internal) or "ext" (external) sample
#' 
#' @return a dataframe that repeats each variant mac times and retains the gene, functional, case, and group status for each variant

make_long = function(counts, leg, case, group) {
  
  # add information to the counts
  temp = counts %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case=case, group=group)
  
  # remove the monomorphic variants 
  # temp2 = temp %>% filter(mac!=0)
  temp2 = temp %>% filter(ac!=0)
  
  # repeat each variant mac times
  # out = data.frame(lapply(temp2, rep, temp2$mac)) %>% select(-count, -mac, -maf)
  out = data.frame(lapply(temp2, rep, temp2$ac)) %>% select(-ac, -af)
  
  return(out)
}

#' merge_cases
#' 
#' @description
#' Function to merge the case datasets for power and t1e calculations for by gene association
#' 
#' @param cases_power case haplotype file used for power calculation
#' @param cases_t1e case haplotype file used for type I error calculation
#' @param leg legend file
#' @param genes_power a vector of gene names used for which power will be calculated
#' 
#' @return the merged case haplotype file containing only the genes being used to calculate t1e and power

merge_cases = function(cases_power, cases_t1e, leg, genes_power) {
  
  # Add row number and gene column to each hap
  hap_power = cases_power %>% mutate(row = leg$row, gene = leg$gene)
  hap_t1e = cases_t1e %>% mutate(row = leg$row, gene = leg$gene)
  
  # Subset haps to the necessary genes
  power_gene = subset(hap_power, gene %in% genes_power) 
  t1e_gene = subset(hap_t1e, !(gene %in% genes_power))
  
  # Merge the two case haps
  hap_out = rbind(power_gene, t1e_gene)
  
  # Order the merged hap file by row number
  hap_out = hap_out[order(hap_out$row),]
  
  # Remove the row number and gene columns
  hap_out = subset(hap_out, select = -c(row, gene))
  
  return(hap_out)
}

#' calc_allele_freqs
#' 
#' @description
#' Function to calculate the ACs/AFs and MACs/MAFs for a given dataset
#' 
#' @param geno a genotype or haplotype matrix denoting the number of reference alleles observed at each individual/haplotype for all variants in the region
#' @param n the number of individuals in geno
#' @param Pop a three letter string denoting the population of geno if applicable (mainly used for reference data)
#' 
#' @return a dataframe denoting the ac, af, mac, and maf for each variant in the region

calc_allele_freqs = function(geno, n, Pop=NULL) {
  
  # counts = data.frame(count = rowSums(geno)) %>%
  #   mutate(mac = ifelse(count>n, 2*n-count, count)) %>%
  #   mutate(maf = mac/(2*n))
  
  # Adelle's way
  counts = data.frame(ac = rowSums(geno)) %>%
    mutate(af = ac/(2*n)) %>%
    mutate(mac = ifelse(ac>n, 2*n-ac, ac)) %>%
    mutate(maf = mac/(2*n))
  
  if(!is.null(Pop)) {
    Pop <- tolower(Pop)
    colnames(counts) <- c(paste0("ac_", Pop), paste0("af_", Pop), paste0("mac_", Pop), paste0("maf_", Pop))
  }

  return(counts)
}

#' est_props
#' 
#' @description
#' Function to estimate the ancestry proportions of a sample using only the common variants
#' 
#' @param counts dataframe containing the ACs and AFs of a sample for each variant in the region
#' @param Pop1 a three letter string denoting the first population in the sample
#' @param Pop2 a three letter string denoting the second population in the sample
#' @param maf a numeric value denoting the minor allele frequency threshold that distinguishes rare variants from common variants
#' 
#' @return returns the proportion estimates outputted by Summix

est_props = function(counts, Pop1, Pop2, maf) {
  
  Pop1 <- tolower(Pop1)
  Pop2 <- tolower(Pop2)
  
  # variants that are common in at least one dataset
  # common <- which(counts$maf > maf | counts[, paste0("maf_", Pop1)] > maf | counts[, paste0("maf_", Pop2)] > maf)
  
  
  # Adelle's way
  # need to filter for both sides of the maf
  common <- which((counts$af > maf & counts$af < 1-maf) | 
                    (counts[, paste0("af_", Pop1)] > maf & counts[, paste0("af_", Pop1)] < 1-maf) |
                    (counts[, paste0("af_", Pop2)] > maf & counts[, paste0("af_", Pop2)] < 1-maf))
  
  
  # Subset counts dataframe to only common variants
  common_df <- counts[common,]
  
  # Use summix to calculate ancestry proportion estimates
  # prop_est <- summix(data = common_df,
  #                    reference=c(paste0("maf_", Pop1), #AFR
  #                                paste0("maf_", Pop2)), #NFE
  #                    observed="maf") #leave out pi.start argument
  
  # Adelle's way
  prop_est <- summix(data = common_df,
                     reference=c(paste0("af_", Pop1), #AFR
                                 paste0("af_", Pop2)), #NFE
                     observed="af", 
                     goodness.of.fit = TRUE, 
                     override_removeSmallRef = TRUE) #show estimates for anc w/ <1% AFs
  
  return(prop_est)
}

#' calc_adjusted_AF
#' 
#' @description
#' Function that uses Summix to update the AFs and ACs of a dataset (primarily the common controls)
#' 
#' @param counts dataframe containing the ACs and AFs for the data to be adjusted as well as the reference data for both populations
#' @param Pop1 a three letter string denoting the first population in the sample
#' @param Pop2 a three letter string denoting the second population in the sample
#' @param case_est the proportion estimates from Summix for the cases
#' @param control_est the proportion estimates from Summix for the controls
#' @param Nref the number of individuals in the reference populations (assumes the number is the same for all reference pops)
#' @param Ncc the number of individuals in the (common) controls
#' 
#' @return a dataframe with the adjusted MACs and MAFs for the controls for each variant in the region

calc_adjusted_AF = function(counts, Pop1, Pop2, case_est, control_est, Nref, Ncc) {
  
  Pop1 <- tolower(Pop1)
  Pop2 <- tolower(Pop2)
  
  # adj_AF <- adjAF(data = counts,
  #                 reference = c(paste0("maf_", Pop1), paste0("maf_", Pop2)),
  #                 observed = "maf",
  #                 pi.target = c(pi_tar1, pi_tar2), 
  #                 pi.observed = c(prop_est[, paste0("maf_", Pop1)], prop_est[, paste0("maf_", Pop2)]),
  #                 adj_method = "average",
  #                 N_reference = c(Nref, Nref),
  #                 N_observed = Ncc,
  #                 filter = TRUE)
  
  # Adelle's way
  adj_AF <- adjAF(data = counts,
                  reference = c(paste0("af_", Pop1), paste0("af_", Pop2)),
                  observed = "af",
                  pi.target = c(case_est[, paste0("af_", Pop1)], case_est[, paste0("af_", Pop2)]), 
                  pi.observed = c(control_est[, paste0("af_", Pop1)], control_est[, paste0("af_", Pop2)]),
                  adj_method = "average",
                  N_reference = c(Nref, Nref),
                  N_observed = Ncc,
                  filter = TRUE) 
  
  # Add adj AF to dataframe
  counts$adj_af <- adj_AF$adjusted.AF$adjustedAF
  
  # Add adj AC to dataframe
  counts$adj_ac <- round(counts$adj_af*(2*Ncc))
  
  # Calculate the adjusted MINOR AF
  counts$adj_maf <- ifelse(counts$adj_af > .5, 1-counts$adj_af, counts$adj_af)
  
  # Calculate the adjusted Minor AC 
  counts$adj_mac <- round(counts$adj_maf*(2*Ncc))
  
  # Return just the adjusted data
  counts_adj <- counts[, c("adj_ac", "adj_af", "adj_mac", "adj_maf")]
  
  # Rename columns so they are same as other data frames
  colnames(counts_adj) <- c("ac", "af", "mac", "maf")
  
  return(counts_adj)
}

# version that doesn't deselect count since the adjusted count files don't have that col
# make_long_adj = function(counts, leg, case, group) {
#   
#   # add information to the counts
#   temp = counts %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case=case, group=group)
#   
#   # remove the monomorphic variants 
#   temp2 = temp %>% filter(mac!=0)
#   
#   # repeat each variant mac times
#   out = data.frame(lapply(temp2, rep, temp2$mac)) %>% select(-mac, -maf)
#   
#   return(out)
# }

### Determine the number of rare fun and syn minor alleles in a dataset FROM THE COUNTS DF
# rare_counts = function(counts, leg.fun, leg.syn, maf){
#   
#   fun.counts = counts[leg.fun$row, ]
#   rare.fun = which(fun.counts$maf <= maf)
#   out.fun = sum(fun.counts[rare.fun, ]$mac)
#   
#   syn.counts = counts[leg.syn$row, ]
#   rare.syn = which(syn.counts$maf <= maf)
#   out.syn = sum(syn.counts[rare.syn, ]$mac)
#   
#   out = c(out.fun, out.syn)
#   
#   return(out)
# }

# Function to calculate allele counts/freqs for all datasets
# calc_allele_freqs_all = function(counts_cases, counts_int, counts_cc, Ncase, Nint, Ncc) {
#   
#   # counts_all = data.frame(count = counts_cases$count + counts_int$count + counts_cc$count) %>% 
#   #   mutate(mac = ifelse(count > (Ncase+Nint+Ncc), 2*(Ncase+Nint+Ncc)-count, count)) %>%
#   #   mutate(maf = mac/(2*(Ncase+Nint+Ncc)))
#   
#   # Adelle's way
#   counts_all = data.frame(ac = counts_cases$ac + counts_int$ac + counts_cc$ac) %>% 
#     mutate(af = ac/(2*(Ncase+Nint+Ncc)))
#   
#   return(counts_all)
# }

# Function to calculate allele counts/freqs for reference datasets
# calc_allele_freqs_ref = function(Pop, hap_ref, Nref) {
#   
#   # counts_ref = data.frame(count = rowSums(hap_ref)) %>%
#   #   mutate(mac = ifelse(count > Nref, 2*Nref-count, count)) %>%
#   #   mutate(maf = mac/(2*Nref))
#   # 
#   # Pop <- tolower(Pop)
#   # colnames(counts_ref) <- c(paste0("count_", Pop), paste0("mac_", Pop), paste0("maf_", Pop))
#   
#   # Adelle's way
#   counts_ref = data.frame(ac = rowSums(hap_ref)) %>%
#     mutate(af = ac/(2*Nref))
# 
#   Pop <- tolower(Pop)
#   colnames(counts_ref) <- c(paste0("ac_", Pop), paste0("af_", Pop))
#   
#   return(counts_ref)
# }

# Function to calculate adjusted MACs and MAFs
# calc_adj_allele_freqs = function(counts, Ncc) {
#   
#   counts = counts %>% mutate(adj_mac2 = ifelse(adj_mac>Ncc, 2*Ncc-adj_mac, adj_mac)) %>%
#     mutate(adj_maf2 = ifelse(adj_maf>0.5, 1-adj_maf, adj_maf))
#   
#   return(counts)
# }

