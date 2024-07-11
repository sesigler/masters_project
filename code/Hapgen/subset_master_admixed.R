################################################################################
# This code subsets legend files from the master legend files for an admixed pop
# It only needs to be run when the admixture proportions of the admixed pop change
################################################################################

# Function to add Pop suffixes to prob column of dfs in pop_probs 
add_suffix <- function(df, suffix) {
  colnames(df)[-1] <- paste0(colnames(df)[-1], ".", suffix)
  return(df)
}

# Pop1 = 'AFR'
# Pop2 = 'NFE'
# admx_pop1 = 80
# admx_pop2 = 20
admx_pop = 'LTX'
Pops = c('IAM', 'NFE', 'EAS', 'AFR')
admx_props = c(47, 44, 5, 4)
# admx_pop = 'AFR_NFE'
# Pops = c('AFR', 'NFE')
# admx_props = c(80, 20)

end = 100 # change back to mysim
start = end-99
set.seed(12345)

# define the admixture proportion
prop_ad = 0.80

library(dplyr)
library(purrr)

dir_in = '/home/math/siglersa/master_legend_files/'
dir_out = paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/subset_master/')

dir_in = 'C:/Users/sagee/Documents/HendricksLab/admixed/master_legs/'

# read in the master legend files
# master.AFR = read.table(paste0(dir_in, 'chr19.block37.AFR.master.legend'), sep='\t')
# master.NFE = read.table(paste0(dir_in, 'chr19.block37.NFE.master.legend'), sep='\t')
# master.IAM = read.table(paste0(dir_in, 'chr19.block37.IAM.master.legend'), sep='\t')
# master.EAS = read.table(paste0(dir_in, 'chr19.block37.EAS.master.legend'), sep='\t')

# Create empty list to store master legend files
master_leg <- setNames(vector("list", length(Pops)), paste0("master.", Pops))

# read in the master legend files
for (i in 1:length(Pops)) {
  master_leg[[i]] <- read.table(paste0(dir_in, 'chr19.block37.', Pops[i], '.master.legend'), sep='\t')
  
  # rename the column names
  colnames(master_leg[[i]]) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
  
  # create an alleles column in the format of a0/a1
  master_leg[[i]]$alleles = paste0(master_leg[[i]]$a0, '/', master_leg[[i]]$a1)
}

# # rename the column names
# colnames(master.AFR) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
# colnames(master.NFE) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
# colnames(master.IAM) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
# colnames(master.EAS) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
# 
# # create an alleles column in the format of a0/a1
# master.AFR$alleles = paste0(master.AFR$a0, '/', master.AFR$a1)
# master.NFE$alleles = paste0(master.NFE$a0, '/', master.NFE$a1)
# master.IAM$alleles = paste0(master.IAM$a0, '/', master.IAM$a1)
# master.EAS$alleles = paste0(master.EAS$a0, '/', master.EAS$a1)

# Create empty list to store probabilities
pop_probs = setNames(vector("list", length(Pops)), paste0(Pops, ".prob"))

for (i in 1:length(Pops)) {
  
  # extract the probabilities at each position for each ancestry
  pop_probs[[i]] = master_leg[[i]] %>% select(position, prob) %>% unique()
}

# # extract the probabilities at each position for each ancestry
# AFR.prob = master.AFR %>% select(position, prob) %>% unique()
# NFE.prob = master.NFE %>% select(position, prob) %>% unique()
# IAM.prob = master.IAM %>% select(position, prob) %>% unique()
# EAS.prob = master.EAS %>% select(position, prob) %>% unique()

# probs = merge(AFR.prob, NFE.prob, by="position", suffixes=c(".AFR", ".NFE"))

# Apply the suffix function to each dataframe in the list
pop_probs2 <- mapply(add_suffix, pop_probs, Pops, SIMPLIFY = FALSE)

# Merge the dataframes in pop_probs2 by the position column
probs <- Reduce(function(x, y) merge(x, y, by = "position", all = TRUE), pop_probs2)

# Create empty list to store singles, dups, and trips for each pop
singles = setNames(vector("list", length(Pops)), paste0("singles.", Pops))
dups = setNames(vector("list", length(Pops)), paste0("dups.", Pops))
trips = setNames(vector("list", length(Pops)), paste0("trips.", Pops))

# subset master legend files according to number of alleles at each position
for (i in seq_along(Pops)) {
  singles[[i]] = master_leg[[i]] %>% filter(prob==1)
  dups[[i]] = master_leg[[i]] %>% filter(prob==0.5)
  trips[[i]] = master_leg[[i]] %>% filter(prob==".")
}
  
# # subset the AFR file according to the number of alleles at each position
# singles.AFR = master.AFR %>% filter(prob==1)
# dups.AFR = master.AFR %>% filter(prob==0.5)
# trips.AFR = master.AFR %>% filter(prob==".")
#   
# # subset the NFE file according to the number of alleles at each position
# singles.NFE = master.NFE %>% filter(prob==1)
# dups.NFE = master.NFE %>% filter(prob==0.5)
# trips.NFE = master.NFE %>% filter(prob==".")
# 
# # subset the AFR file according to the number of alleles at each position
# singles.IAM = master.IAM %>% filter(prob==1)
# dups.IAM = master.IAM %>% filter(prob==0.5)
# trips.IAM = master.IAM %>% filter(prob==".")
# 
# # subset the NFE file according to the number of alleles at each position
# singles.EAS = master.EAS %>% filter(prob==1)
# dups.EAS = master.EAS %>% filter(prob==0.5)
# trips.EAS = master.EAS %>% filter(prob==".")

# Create empty list to store single positions seen in all populations
singles2 = setNames(vector("list", length(Pops)), paste0("singles.", Pops, "2"))

### Problem here: length of singles differ between pops after filtering
# subset the single positions seen in all populations
for (i in seq_along(singles)) {
  
  # Concatenate the position columns of all dfs in the list except the current df
  other_positions <- unlist(lapply(singles[-i], function(df) df$position))
  
  # Get unique positions from the concatenated positions
  unique_positions <- unique(other_positions)
  
  # Filter the current df to the positions seen in all pops
  singles2[[i]] = singles[[i]] %>% filter(position %in% unique_positions)

}
  
# subset the single positions seen in both populations
# singles.AFR2 = singles.AFR %>% filter(position %in% singles.NFE$position)
# singles.NFE2 = singles.NFE %>% filter(position %in% singles.AFR$position)

singles.AFR2 = singles.AFR %>% filter(position %in% c(singles.NFE$position, singles.IAM$position, singles.EAS$position))
singles.NFE2 = singles.NFE %>% filter(position %in% c(singles.AFR$position, singles.IAM$position, singles.EAS$position))
singles.IAM2 = singles.IAM %>% filter(position %in% c(singles.NFE$position, singles.AFR$position, singles.EAS$position))
singles.EAS2 = singles.EAS %>% filter(position %in% c(singles.NFE$position, singles.IAM$position, singles.AFR$position))

# keep the single positions with same alt allele in all pops
common_alleles <- singles2[-1] %>% map(~ .[["alleles"]]) %>% reduce(intersect)
singles.keep = singles2[[1]] %>% filter(alleles==singles2$singles.NFE2$alleles)
  
# keep the single positions with the same alternate allele in both populations
# singles.keep = singles.AFR2 %>% filter(alleles==singles.NFE2$alleles)
singles.keep.OG = singles.AFR2 %>% filter(alleles==singles.NFE2$alleles)
  
# subset the single positions in just one of the populations
singles.NA.NFE = singles.NFE %>% filter(!(position %in% singles.AFR$position))
singles.NA.AFR = singles.AFR %>% filter(!(position %in% singles.NFE$position))
  
# combine the single positions together
singles.keep2 = rbind(singles.keep, singles.NA.NFE, singles.NA.AFR)
  
# subset the duplicate positions in both populations
dups.po = probs %>% filter((prob.AFR==0.5 & prob.NFE==0.5) | (prob.AFR==0.5 & prob.NFE==".") | (prob.AFR=="." & prob.NFE==0.5)) %>% select(position)
dups.AFR2 = dups.AFR %>% filter(position %in% dups.po$position)
dups.NFE2 = dups.NFE %>% filter(position %in% dups.po$position) %>% filter(!(position %in% dups.AFR2$position))
dups.keep = rbind(dups.AFR2, dups.NFE2)

# subset the triplicate positions in both populations
trips.po = probs %>% filter(prob.AFR=="." & prob.NFE==".") %>% select(position)
trips = trips.AFR %>% filter(position %in% trips.po$position)
    
for (i in start:end){
    
  # create a table of transitions/transversions for the triplicates
  trips.po1  =  as.data.frame(matrix(NA, nrow=nrow(trips.po), ncol=3))
  colnames(trips.po1) = c('position', 'draw', 'TiTv')
  trips.po1$position = trips.po$position
  trips.po1$draw = runif(nrow(trips.po1))
  trips.po1$TiTv[which(trips.po1$draw < 0.7396)] = 'transition'
  trips.po1$TiTv[which(trips.po1$draw >= 0.7396)] = 'transversion'
  #head(trips.po1)
    
  # subset the transitions
  ti = trips.po1 %>% filter(TiTv  == 'transition')
  ann.ti = trips %>% filter(position %in% ti$position, alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C')
    
  # remove all of the transitions
  ann.tv = trips %>% filter(!(position %in% ti$position), !(alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C'))
    
  # randomly pick an allele from the transversions
  ann.tv2 = ann.tv %>% group_by(position) %>% sample_n(size=1)
    
  # merge transitions and transversions
  trips2 = union(ann.ti, ann.tv2)
    
  # randomly pick an allele from the duplicates
  dups2 = dups.keep %>% group_by(position) %>% sample_n(size=1)
    
  # subset the single positions with different alternate alleles in the populations
  singles.dif = singles.AFR2 %>% filter(alleles!=singles.NFE2$alleles)
    
  # pick an allele for the single positions with different alternate alleles based on the admixed proportion
  singles.dif.AFR = singles.dif %>% sample_n(size=round(prop_ad*nrow(singles.dif)))
  singles.dif.NFE = singles.NFE2 %>% filter(alleles!=singles.AFR2$alleles) %>% filter(!(position %in% singles.dif.AFR$position))
    
  # combine the single positions together
  singles2 = rbind(singles.keep2, singles.dif.AFR, singles.dif.NFE)
    
  # merge all positions
  master2 = union(singles2, union(dups2, trips2)) %>% arrange(position) %>% 
    select(id, position, a0, a1, AC, prob, exonic, gene, fun)

  master2$fun = ifelse(master2$fun=="synonymous SNV", "syn", "fun")
  master2$exonic[grepl("exonic", master2$exonic)] = "exonic"
  master2$gene[grepl("ZNF333", master2$gene)] = "ZNF333"
    
  # write output legend file
  out.name = paste0(dir_out, 'chr19.block37.', admx_pop, '.sim', i, '.copy.legend')
  write.table(master2, out.name, row.names=F, col.names=T, quote=F, sep='\t')

  print(i)
  #print(nrow(master2))
}