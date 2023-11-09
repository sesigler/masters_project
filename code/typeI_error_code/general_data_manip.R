############################################################################## 
# This file contains the functions necessary to do general data manipulation
# on the files necessary for performing type I error and power calculations for 
# several different rare variant association tests
#
# Variable legend:
# hap: the haplotype file
# counts: dataframe containing the allele counts, MACs, and MAFs of a genotype file
# leg: the legend file
# case: a string indicating if the counts are "cases" or "controls"
# group: a string indicating if the counts are from an "int" (internal) or 
#        "ext" (external) sample
# geno: the genotype file
# n: number of individuals in the dataset/number of columns in the genotype file
# N<DATASET>: number of individuals in DATASET
# Pop: a string referring to the 3 letter ancestry population (e.g. "AFR", "NFE")
# hap_ref: haplotype file for the reference dataset
# maf: the minor allele frequency threshold 
# prop_est: ancestry proportion estimates from the summix function
# pi_tar: the target proportion of a given ancestry
# cases_power: hap file for cases used for power calculation
# cases_t1e: hap file for cases used for type I error calculation
# power_genes: vector of gene names used for power calculation
##############################################################################

# function to convert the haplotypes into a genotype matrix
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

# function to create a dataframe with a line for each variant instead of just counts 
# (necessary for ProxECAT v2)
# version that doesn't filter out common variants just yet
make_long = function(counts, leg, case, group) {
  
  # add information to the counts
  temp = counts %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case=case, group=group)
  
  # remove the monomorphic variants 
  temp2 = temp %>% filter(mac!=0)
  
  # repeat each variant mac times
  out = data.frame(lapply(temp2, rep, temp2$mac)) %>% select(-count, -mac, -maf)
  
  return(out)
}

# version that doesn't deselect count since the adjusted count files don't have that col
make_long_adj = function(counts, leg, case, group) {
  
  # add information to the counts
  temp = counts %>% mutate(id=leg$id, gene=leg$gene, fun=leg$fun, case=case, group=group)
  
  # remove the monomorphic variants 
  temp2 = temp %>% filter(mac!=0)
  
  # repeat each variant mac times
  out = data.frame(lapply(temp2, rep, temp2$mac)) %>% select(-mac, -maf)
  
  return(out)
}

# Merge the case datasets for power and t1e calculations
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

# Function to calculate the allele counts/frequencies
calc_allele_freqs = function(geno, n) {
  
  counts = data.frame(count = rowSums(geno)) %>% 
    mutate(mac = ifelse(count>n, 2*n-count, count)) %>%
    mutate(maf = mac/(2*n))
  
  return(counts)
}

# Function to calculate allele counts/freqs for all datasets
calc_allele_freqs_all = function(counts_cases, counts_int, counts_cc, Ncase, Nint, Ncc) {
  
  counts_all = data.frame(count = counts_cases$count + counts_int$count + counts_cc$count) %>% 
    mutate(mac = ifelse(count > (Ncase+Nint+Ncc), 2*(Ncase+Nint+Ncc)-count, count)) %>%
    mutate(maf = mac/(2*(Ncase+Nint+Ncc)))
  
  return(counts_all)
}

# Function to calculate allele counts/freqs for reference datasets
calc_allele_freqs_ref = function(Pop, hap_ref, Nref) {
  
  counts_ref = data.frame(count = rowSums(hap_ref)) %>% 
    mutate(mac = ifelse(count > Nref, 2*Nref-count, count)) %>%
    mutate(maf = mac/(2*Nref))
  
  Pop <- tolower(Pop)
  colnames(counts_ref) <- c(paste0("count_", Pop), paste0("mac_", Pop), paste0("maf_", Pop))
  
  return(counts_ref)
}

# Function to calculate adjusted MACs and MAFs
calc_adj_allele_freqs = function(counts, Ncc) {
  
  counts = counts %>% mutate(adj_mac2 = ifelse(adj_mac>Ncc, 2*Ncc-adj_mac, adj_mac)) %>%
    mutate(adj_maf2 = ifelse(adj_maf>0.5, 1-adj_maf, adj_maf))
  
  return(counts)
}

# Estimate ancestry proportions using only common variants
est_props = function(counts, Pop1, Pop2, maf) {
  
  Pop1 <- tolower(Pop1)
  Pop2 <- tolower(Pop2)
  
  # variants that are common in at least one dataset
  common <- which(counts$maf > maf | counts[, 6] > maf | counts[, 9] > maf) #6=maf_afr, 9=maf_nfe
  
  # Subset counts dataframe to only common variants
  common_df <- counts[common,]
  
  # Use summix to calculate ancestry proportion estimates
  prop_est <- summix(data = common_df,
                     reference=c(paste0("maf_", Pop1), #AFR
                                 paste0("maf_", Pop2)), #NFE
                     observed="maf") #leave out pi.start argument
  
  return(prop_est)
}

# Use summix to update AFs of common controls dataset
calc_adjusted_AF = function(counts, Pop_ref, prop_est, pi_tar1, pi_tar2, Ncc) {
  
  Pop_ref <- tolower(Pop_ref)
  
  adj_AF <- adjAF(data = counts,
                    reference = c(paste0("maf_", Pop_ref)),
                    observed = "maf",
                    pi.target = c(pi_tar1, pi_tar2), #last one is AFR proportion
                    pi.observed = c(prop_est[, 6], prop_est[, 5])) #6=maf_nfe, 5=maf_afr
  
  # Add adj AF to data frame
  counts$adj_maf <- adj_AF$adjusted.AF$adjustedAF
  
  # Set AFs < 0 = 0 
  counts$adj_maf[counts$adj_maf < 0] <- 0
  
  # Add adj ACs
  counts$adj_mac <- round(counts$adj_maf*(2*Ncc))
  
  # Re-check that ACs and AFs are the minor ones (may 68be higher after adjustment)
  counts_minor = calc_adj_allele_freqs(counts, Ncc)
  
  # Create data frame with only the 2 adj MAC and AF columns
  counts_adj <- counts_minor[, c("adj_mac2", "adj_maf2")]
  
  # Rename columns so they are same as other data frames
  colnames(counts_adj) <- c("mac", "maf")
  
  return(counts_adj)
}

