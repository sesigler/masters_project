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


#' flip_file
#' 
#' @description
#' Function the flip values for a specified file at variants with an AF >= 1-maf
#' 
#' @param file_to_flip the file to be updated
#' @param flip the rows of the legend file specifying which variants need to be flipped
#' @param N number of individuals, only for flipping adjusted common controls
#' @param file_type string specifying what file is to be flipped
#' 
#' @return the specified file with the relevant values flipped at the variants in flip

flip_file = function(file_to_flip, flip, N=NULL, file_type) {
  
  # Make copy of original file
  file2 = file_to_flip
  
  if (file_type == "leg") {
    
    # Flip ref allele in file2 with alt allele in file
    file2$a0[flip$row] <- file_to_flip$a1[flip$row]
    
    # Flip alt allele in leg2 with ref allele in leg
    file2$a1[flip$row] <- file_to_flip$a0[flip$row]
    
    return(file2)
    
  } else if (file_type == "geno") {
    
    # Flip the alternate allele counts at the relevant variants
    file2[flip$row,] <- 2-file_to_flip[flip$row,]
    
    return(file2)
    
  } else if (file_type == "count") {
    
    # Flip the AFs at the variants that need to be flipped
    file2$af[flip$row] <- 1-file_to_flip$af[flip$row]
    
    # Update the ACs of the variants that were flipped
    file2$ac[flip$row] <- round(file2$af[flip$row]*(2*N))
    
    return(file2)
    
  } else {
    stop("ERROR: 'file_type' must be a string of either 'leg', 'geno', or 'count'")
  }
  
}

#' flip_data
#' 
#' @description
#' Function to flip the relevant datasets at variants where AF >= 1-maf for each possible scenario
#' 
#' @param leg legend file
#' @param flip the rows of the legend file specifying which variants need to be flipped
#' @param geno.case genotype matrix for the cases
#' @param geno.int genotype matric for the internal controls, default value is NULL
#' @param geno.cc genotype matric for the common controls, default value is NULL
#' @param count.cc.adj adjusted AC and AF dataframe for the commcon controls, default value is NULL
#' @param Ncase number of individuals in the cases
#' @param Nint number of individuals in the internal controls, default value is NULL
#' @param Ncc number of individuals in the common controls, default value is NULL
#' @param cntrl string specifying the type(s) of controls used, options are int, ext, and all
#' @param adj boolean value specifying if the common controls used are adjusted or unadjusted, default value is false
#' 
#' @return a list of all the relevant files that could be flipped, returns the unaltered files if no variants needed to be flipped

flip_data = function(leg, flip, geno.case, geno.int=NULL, geno.cc=NULL, count.cc.adj=NULL, Ncase, Nint=NULL, Ncc=NULL, cntrl, adj=FALSE) {
  
  if (nrow(flip) != 0) {
    
    if (cntrl == "int") {
      
      # Create new leg file
      leg2 = flip_file(leg, flip, N=NULL, file_type="leg") 
      
      # Update geno files
      geno_cases2 = flip_file(geno_cases, flip, N=NULL, file_type="geno")
      geno_int2 = flip_file(geno_int, flip, N=NULL, file_type="geno")
      
      # Recalculate ac/af 
      count_cases2 = calc_allele_freqs(geno_cases2, Ncase, Pop=NULL)
      count_int2 = calc_allele_freqs(geno_int2, Nint, Pop=NULL)
      
      # Return all changed files
      return(list(leg2, geno_cases2, geno_int2, count_cases2, count_int2))
      
    } else if (cntrl == "ext" & !adj) {
      
      # Create new leg file
      leg2 = flip_file(leg, flip, N=NULL, file_type="leg")
      
      # Update geno files
      geno_cases2 = flip_file(geno_cases, flip, N=NULL, file_type="geno")
      geno_cc2 = flip_file(geno_cc, flip, N=NULL, file_type="geno")
      
      # Recalculate ac/af 
      count_cases2 = calc_allele_freqs(geno_cases2, Ncase, Pop=NULL)
      count_cc2 = calc_allele_freqs(geno_cc2, Ncc, Pop=NULL)
      
      # Return all changed files
      return(list(leg2, geno_cases2, geno_cc2, count_cases2, count_cc2))
      
    } else if (cntrl == "ext" & adj) {
      
      # Create new leg file
      leg2 = flip_file(leg, flip, N=NULL, file_type="leg")
      
      # Update geno files
      geno_cases2 = flip_file(geno_cases, flip, N=NULL, file_type="geno")
      
      # Recalculate ac/af 
      count_cases2 = calc_allele_freqs(geno_cases2, Ncase, Pop=NULL)
      count_cc_adj2 = flip_file(count_cc_adj, flip, N=Ncc, file_type="count")
      
      # Return all changed files
      return(list(leg2, geno_cases2, count_cases2, count_cc_adj2))
      
    } else if (cntrl == "all" & !adj) {
      
      # Create new leg file
      leg2 = flip_file(leg, flip, N=NULL, file_type="leg")
      
      # Update geno files
      geno_cases2 = flip_file(geno_cases, flip, N=NULL, file_type="geno")
      geno_int2 = flip_file(geno_int, flip, N=NULL, file_type="geno")
      geno_cc2 = flip_file(geno_cc, flip, N=NULL, file_type="geno")
      
      # Recalculate ac/af 
      count_cases2 = calc_allele_freqs(geno_cases2, Ncase, Pop=NULL)
      count_int2 = calc_allele_freqs(geno_int2, Nint, Pop=NULL)
      count_cc2 = calc_allele_freqs(geno_cc2, Ncc, Pop=NULL)
      
      # Return all changed files
      return(list(leg2, geno_cases2, geno_int2, geno_cc2, count_cases2, count_int2, count_cc2))
      
    } else if (cntrl == "all" & adj) {
      
      # Create new leg file
      leg2 = flip_file(leg, flip, N=NULL, file_type="leg")
      
      # Update geno files
      geno_cases2 = flip_file(geno_cases, flip, N=NULL, file_type="geno")
      geno_int2 = flip_file(geno_int, flip, N=NULL, file_type="geno")
      
      # Recalculate ac/af 
      count_cases2 = calc_allele_freqs(geno_cases2, Ncase, Pop=NULL)
      count_int2 = calc_allele_freqs(geno_int2, Nint, Pop=NULL)
      count_cc_adj2 = flip_file(count_cc_adj, flip, N=Ncc, file_type="count")
      
      
      # Return all changed files
      return(list(leg2, geno_cases2, geno_int2, count_cases2, count_int2, count_cc_adj2))
      
    }
    
  } else if (cntrl == "int") {
    
    # Return all relevant files
    return(list(leg, geno_cases, geno_int, count_cases, count_int))
    
  } else if (cntrl == "ext" & !adj) {
    
    return(list(leg, geno_cases, geno_cc, count_cases, count_cc))
    
  } else if (cntrl == "ext" & adj) {
    
    return(list(leg, geno_cases, count_cases, count_cc))
    
  } else if (cntrl == "all" & !adj) {
    
    return(list(leg, geno_cases, geno_int, geno_cc, count_cases, count_int, count_cc))
    
  } else if (cntrl == "all" & adj) {
    
    return(list(leg, geno_cases, geno_int, count_cases, count_int, count_cc))
  }
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

