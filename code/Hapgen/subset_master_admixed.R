################################################################################
# This code subsets rows from the master legend files to create legend files
# for an admixed pop
# It only needs to be run when the admixture proportions of the admixed pop change
################################################################################

library(dplyr)
library(purrr)

# Define parameters
admx_pop <- 'LTX'
Pops <- c('IAM', 'NFE', 'EAS', 'AFR') 
admx_props <-  c(IAM = 47, NFE = 44, EAS = 5, AFR = 4)

end <- 100 # change back to mysim
start <- end-99
set.seed(12345)

# Directory paths
dir_in <- '/home/math/siglersa/master_legend_files/'
dir_out <- paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/subset_master/')

# dir_in <- 'C:/Users/sagee/Documents/HendricksLab/admixed/master_legs/'


# Create empty lists to store master legend files, probabilities, singles, dups, and trips for each pop
master_leg <- setNames(vector("list", length(Pops)), paste0("master.", Pops))
pop_probs <- setNames(vector("list", length(Pops)), paste0(Pops, ".prob"))
singles <- setNames(vector("list", length(Pops)), paste0("singles.", Pops))
dups <- setNames(vector("list", length(Pops)), paste0("dups.", Pops))
trips <- setNames(vector("list", length(Pops)), paste0("trips.", Pops))

# Loop through each pop
for (pop in 1:length(Pops)) {
  
  # Read in legend file
  master_leg[[pop]] <- read.table(paste0(dir_in, 'chr19.block37.', Pops[pop], '.master.legend'), sep='\t')
  
  # rename the column names
  colnames(master_leg[[pop]]) <- c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
  
  # create an alleles column in the format of a0/a1
  master_leg[[pop]]$alleles <- paste0(master_leg[[pop]]$a0, '/', master_leg[[pop]]$a1)
  
  # extract the probabilities at each position for each ancestry
  new_name = paste0("prob.", Pops[pop])
  pop_probs[[pop]] <- master_leg[[pop]] %>% select(position, prob) %>% rename(!!new_name := prob) %>% unique()
  
  # subset master legend files according to number of alleles at each position
  singles[[pop]] <- master_leg[[pop]] %>% filter(prob==1)
  dups[[pop]] <- master_leg[[pop]] %>% filter(prob==0.5)
  trips[[pop]] <- master_leg[[pop]] %>% filter(prob==".")
}

# Merge the dataframes in pop_probs by the position column
probs <- Reduce(function(x, y) merge(x, y, by = "position", all = TRUE), pop_probs)

# Calculate the number of 1s, 0.5s, and '.'s in each row, excluding the position column
num_ones <- rowSums(probs[, -1] == 1)
num_half <- rowSums(probs[, -1] == 0.5)
num_dot <- rowSums(probs[, -1] == '.')

### DIFFERENT SINGLE SCENARIOS
# 1. Prob = 1 in all pops and they all have same alt allele (choose one pop to subset positions from its legend)
# 2. Prob = 1 in all pops but they have diff alt alleles (use admixture props to choose alt allele)
# 3. Prob = 1 in only one pop (keep that pop's alt allele)
# 4. Prob = 1 in more than one pop but not every pop (again use admixture props, but only select from pops whose prob = 1)

### Singles scenarios 1 and 2

# Filter probs df to rows where all probability columns have value of 1
probs_all1 <- probs[num_ones == length(Pops), ]

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

# Subset first pop leg file to positions in singles_same_alt
# Note: ref and alt allele same, so doesn't matter which pop's leg file you subset the rows from
singles.same <- singles[[1]] %>% filter(position %in% singles_same_alt)

# Note: rows for singles scenario 2 are subset below in main for loop

### Singles scenario 3

# Filter probs df to rows where only one pop has prob = 1
probs_one1 = probs[num_ones == 1, ]

# Create a vector to store row from leg file where only one pop has prob = 1
singles.one1 <- c()

# Loop through each position in probs_one1 to get position where only that pop has prob = 1
for (pos in probs_one1$position) {

  # Obtain the col name which has prob = 1 at current position and remove "prob." from the name
  current_pop <- gsub("prob\\.", "", names(probs_one1[probs_one1$position == pos, ])[which(probs_one1[probs_one1$position == pos, ] == 1)])
  
  # Subset the row from the current pop leg file to the current position
  leg_row <- filter(singles[[paste0('singles.', current_pop)]], position==pos)

  # Add the row from the current pop's legend file to the singles.one1 vector
  singles.one1 <- rbind(singles.one1, leg_row)

}

### Singles scenario 4 (rows subset below in main for loop)
# Note: there are other ways you can do this, but this was the easiest

# Filter probs df to positions where more than one pop has prob = 1 but not every pop has prob = 1
probs_not_all1 <- probs[num_ones > 1 & num_ones < length(Pops), ]


### DUPLICATE SCENARIO (rows subset below in main for loop)
# Similar to scenario 4 from singles

# Filter probs df to positions where at least one pop has prob = 0.5 but no pop has prob = 1
probs_dups_no1 <- probs[num_ones == 0 & num_half > 0, ]


### TRIPLICATES 

# Filter probs df to triplicate positions across all populations
trips.po <- probs[num_dot == length(Pops), ]

# Create vector to store same triplicate positions across all pops
# All positions and alleles are same across all pops so just subset from first pop
trips.pop <- trips[[1]] %>% filter(position %in% trips.po$position)

    
# Main for loop to add randomness to admixed legend file for each simulation rep
for (i in start:end){
    
  ### TRIPLICATES 
  
  # create a table of transitions/transversions for the triplicates
  trips.po1  =  as.data.frame(matrix(NA, nrow=nrow(trips.po), ncol=3))
  colnames(trips.po1) = c('position', 'draw', 'TiTv')
  trips.po1$position = trips.po$position
  trips.po1$draw = runif(nrow(trips.po1))
  trips.po1$TiTv[which(trips.po1$draw < 0.7396)] = 'transition'
  trips.po1$TiTv[which(trips.po1$draw >= 0.7396)] = 'transversion'
    
  # subset the transitions
  ti = trips.po1 %>% filter(TiTv  == 'transition')
  
  # subsetting from the first pop legend file since we filtered to the common 
  # triplicate positions so it doesn't matter which legend file we subset from
  ann.ti = trips.pop %>% filter(position %in% ti$position, alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C')
    
  # remove all of the transitions
  ann.tv = trips.pop %>% filter(!(position %in% ti$position), !(alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C'))
    
  # randomly pick an allele from the transversions
  ann.tv2 = ann.tv %>% group_by(position) %>% sample_n(size=1)
    
  # merge transitions and transversions
  trips2 = union(ann.ti, ann.tv2)
    
  ### DUPLICATES 
  
  # Create vector to store duplicate rows for each pop 
  dups_keep <- c()
  
  # Loop through each row and choose a population for each position based on the admixture proportions
  for (j in 1:nrow(probs_dups_no1)) {
    
    # Subset dataframe to current row
    row <- probs_dups_no1[j, ]
    
    # select only the pops with prob = 0.5 at current row
    available_pops <- Pops[which(row[-1] == 0.5)] 
    
    # Use admixture props to choose which population's allele info to use for current position
    # Note: sample will first normalize the probabilities if they don't sum to 1
    selected_pop <- sample(available_pops, 1, prob = admx_props[available_pops]/100) 
    
    # Subset the selected pop's leg file to the current position
    pop_rows <- filter(dups[[paste0('dups.', selected_pop)]], position==row$position)
    
    # Add rows to vector
    dups_keep <- rbind(dups_keep, pop_rows)
    
  }
  
  # randomly pick an allele from the duplicates
  dups2 = dups_keep %>% group_by(position) %>% sample_n(size=1)
  
  ### SINGLES
  
  ### Singles scenario 2
  
  # Create a vector to store the rows of the legend files chosen for each position in singles_diff_alt
  singles.diff <- c()
  
  # Loop through each position
  for (pos in singles_diff_alt) {
    
    # pick an allele for the single positions with diff alt alleles based on the admixed proportions
    selected_pop <- sample(Pops, 1, prob = admx_props/100)
    
    # Subset the selected pop's leg file to the current position
    pop_row <- filter(singles[[paste0('singles.', selected_pop)]], position==pos)
    
    # Add the row from the selected pop's legend file to the singles diff vector
    singles.diff <- rbind(singles.diff, pop_row)
  }
  
  ### Singles scenario 4
  
  # Create a list to store the rows of the legend files chosen for each position in singles_not_all1
  singles.some1 <- c()
  
  # Loop through each row and choose a population for each position based on the admixture proportions
  for (k in 1:nrow(probs_not_all1)) {
    
    # Subset dataframe to current row
    row <- probs_not_all1[k, ]
    
    # select only the pops with prob = 1 at current row
    available_pops <- Pops[which(row[-1] == 1)] 
    
    # Use admixture props to choose which population's allele info to use for current position
    selected_pop <- sample(available_pops, 1, prob = admx_props[available_pops]/100) 
    
    # Subset the selected pop's leg file to the current position
    pop_rows <- filter(singles[[paste0('singles.', selected_pop)]], position==row$position)
    
    # Add row to vector
    singles.some1 <- rbind(singles.some1, pop_rows)
  }
    
  # combine the single positions together
  singles2 = rbind(singles.same, singles.diff, singles.one1, singles.some1)
    
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