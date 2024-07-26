################################################################################
# This code subsets rows from the master legend files to create legend files
# for an admixed pop
# It only needs to be run when the admixture proportions of the admixed pop change
################################################################################

# Pop1 = 'AFR'
# Pop2 = 'NFE'
# admx_pop1 = 80
# admx_pop2 = 20
admx_pop <- 'LTX'
Pops <- c('IAM', 'NFE', 'EAS', 'AFR') 
admx_props <-  c(IAM = 47, NFE = 44, EAS = 5, AFR = 4)
# admx_pop = 'AFR_NFE'
# Pops = c('AFR', 'NFE')
# admx_props = c(80, 20)

end <- 100 # change back to mysim
start <- end-99
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

# Create empty lists to store master legend files, probabilities, singles, dups, and trips for each pop
master_leg <- setNames(vector("list", length(Pops)), paste0("master.", Pops))
pop_probs = setNames(vector("list", length(Pops)), paste0(Pops, ".prob"))
singles = setNames(vector("list", length(Pops)), paste0("singles.", Pops))
dups = setNames(vector("list", length(Pops)), paste0("dups.", Pops))
trips = setNames(vector("list", length(Pops)), paste0("trips.", Pops))

# read in the master legend files
for (i in 1:length(Pops)) {
  master_leg[[i]] <- read.table(paste0(dir_in, 'chr19.block37.', Pops[i], '.master.legend'), sep='\t')
  
  # rename the column names
  colnames(master_leg[[i]]) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
  
  # create an alleles column in the format of a0/a1
  master_leg[[i]]$alleles = paste0(master_leg[[i]]$a0, '/', master_leg[[i]]$a1)
  
  # extract the probabilities at each position for each ancestry
  new_name = paste0("prob.", Pops[i])
  pop_probs[[i]] = master_leg[[i]] %>% select(position, prob) %>% rename(!!new_name := prob) %>% unique()
  
  # subset master legend files according to number of alleles at each position
  singles[[i]] = master_leg[[i]] %>% filter(prob==1)
  dups[[i]] = master_leg[[i]] %>% filter(prob==0.5)
  trips[[i]] = master_leg[[i]] %>% filter(prob==".")
}

# Merge the dataframes in pop_probs by the position column
probs <- Reduce(function(x, y) merge(x, y, by = "position", all = TRUE), pop_probs)

### DIFFERENT SINGLE SCENARIOS
# 1. Prob = 1 in all pops and they all have same alt allele (choose one pop to subset positions from its legend)
# 2. Prob = 1 in all pops but they have diff alt alleles (use admixture props to choose alt allele)
# 3. Prob = 1 in only one pop (keep that pop's alt allele)
# 4. Prob = 1 in more than one pop but not every pop (again use admixture props, but only select from pops whose prob = 1)

### Singles scenarios 1 and 2

# Filter probs df to rows where all probability columns have value of 1
probs_all1 = probs %>% filter(if_all(-position, ~ . == 1))

# Create vectors to store positions
singles_same_alt <- c() # positions where prob=1 across all pops AND alt (and ref) alleles are the same  
singles_diff_alt <- c() # positions where prob=1 across all pops but alt alleles differ (same ref allele)

# Loop through all positions where prob=1 across all pops
for (pos in probs_all1$position) {
  
  # Get the ref/alt alleles for each pop at the current position
  alleles <- sapply(master_leg, function(df) df$alleles[df$position == pos])
  
  # Check if the ref/alt alleles are the same
  if (length(unique(alleles)) == 1) {
    
    # If ref/alt alleles are the same, add the position to singles_same_alt
    singles_same_alt <- c(singles_same_alt, pos)
  }
  else {
    # If the alt alleles differ (assuming ref alleles are the same-I checked, they are),
    # add the position to singles_diff_alt
    singles_diff_alt <- c(singles_diff_alt, pos)
  }
}


### Singles scenario 3

# Filter probs df to rows where only one pop has prob = 1
probs_one1 = probs %>% rowwise() %>% filter(sum(c_across(-position) == 1) == 1) %>% ungroup

# Create a list to store positions where only one pop has prob = 1 for each pop
singles_one1_pop <- setNames(vector("list", length(Pops)), paste0(Pops))

# Loop through each prob.pop col to get position where only that pop has prob = 1
for (pop in Pops) {
  singles_one1_pop[[pop]] <- probs_one1 %>% filter(!!sym(paste0('prob.', pop)) == 1) %>% pull(position)
}

### Singles scenario 4
# Note: there are other ways you can do this, but this was the easiest

# Calculate the number of 1s in each row, excluding the position column
num_ones <- rowSums(probs[, pop_cols] == 1)

# Filter probs df to positions where more than one pop has prob = 1 but not every pop has prob = 1
probs_not_all1 <- probs[num_ones > 1 & num_ones < length(pop_cols), ]

# Create a list to store positions for each population
singles_not_all1 <- setNames(vector("list", length(Pops)), paste0(Pops))

# Loop through each row and choose a population for each position based on the admixture proportions
for (i in 1:nrow(probs_not_all1)) {
  
  # Subset dataframe to current row
  row <- probs_not_all1[i, ]
  available_pops <- Pops[which(row[-1] == 1)] # select only the pops with prob = 1 at current row
  
  # Use admixture props to choose which population's allele info to use for current position
  # Note: sample will first normalize the probabilities if they don't sum to 1
  selected_pop <- sample(available_pops, 1, prob = admx_props[available_pops]/100) 

  # store position in the pop's vector
  singles_not_all1[[selected_pop]] <- c(singles_not_all1[[selected_pop]], row$position) 
}



# subset the single positions seen in both populations
# singles.AFR2 = singles.AFR %>% filter(position %in% singles.NFE$position)
# singles.NFE2 = singles.NFE %>% filter(position %in% singles.AFR$position)
  
# keep the single positions with the same alternate allele in both populations
# singles.keep = singles.AFR2 %>% filter(alleles==singles.NFE2$alleles)
  
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