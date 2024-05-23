# This code is broken up into 3 steps
# Step 1: Extract the IAM (hg38) coordinates and format them to be input into LiftOver 
#         to convert them to hg19 coordinates
# Step 2: Update the IAM legend file to have the hg19 coordinates
# Step 3: Remove positions from IAM leg and hap file not seen in the original 
#         legend file (using NFE_Block37_CDS_ref_added.legend from 
#         /storage/math/projects/compinfo/simulations/input/NFE_blocks since the
#         coding bp do not vary by ancestry).
#         Then format IAM hap and leg files so they include all the positions seen
#         in original legend file (add in filler data at the positions not seen in
#         the IAM data).


# load libraries
library(dplyr)
library(tidyr)
library(data.table)

dir = 'C:/Users/sagee/Documents/HendricksLab/IAM_hap_data/'

### Step 1: Write out positions of IAM legend file in LiftOver format for conversion to hg19 coordinates

# Read in IAM legend file
IAM_leg = read.table(paste0(dir, 'hgdp_1kg_phased_haps_v2_block37_IAM.legend'), header=T, sep='')

# Create column of base pair (bp) positions in correct format for LiftOver
IAM_leg = IAM_leg %>% mutate(pos2 = paste0("chr19:", position, "-", position))

# extract formatted position column
out <- IAM_leg %>% select(pos2)

# save positions as txt file
write.table(out, paste0(dir, "IAM_legend_formatted_positons.txt"), quote=F, row.names=F, col.names=F)


### Step 2: Update IAM leg to have hg19 coordinates

# Read in IAM positions (hg38) converted to hg19 coordinates
iam_hg19 <- read.table(paste0(dir, 'IAM_leg_hg38_to_hg19_LiftOver_positions.bed'), header = F)

# Separate LiftOver positions into separate columns  
hg19_pos = iam_hg19 %>% separate(col="V1", into=c("chr", "pos1", "pos2"), sep = ":|-")

# Create new IAM leg with hg19 positions instead of hg38 coordinates
iam_leg_hg19 <- data.frame(position = hg19_pos$pos2, a0 = IAM_leg$a0, a1 = IAM_leg$a1)

# Make id column and make it first column
iam_leg_hg19 = iam_leg_hg19 %>% mutate(id = paste0("chr19:", position, "_", a0, "_", a1), 
                                       row = 1:nrow(iam_leg_hg19)) %>%
  select(id, everything())


### Step 3: Add in missing coding bp and remove any bp not seen in original legend file 
###         to both IAM legend and hap files

# Read in RAREsim legend file (chr19 block 37 positions are same across ancestry)
leg = read.table(paste0(dir, 'NFE_Block37_CDS_ref_added.legend'), header=T, sep='\t')
leg$row <- 1:nrow(leg)

# Read in IAM hap file
iam_hap = fread(paste0(dir, 'hgdp_1kg_phased_haps_v2_block37_IAM.hap'))
iam_hap = as.data.frame(iam_hap)

# Check for duplicate positions in IAM legend file, then remove them
summary(duplicated(iam_leg_hg19$position)) # 131
iam_leg2 <- iam_leg_hg19[-c(which(duplicated(iam_leg_hg19$position))), ]
iam_leg2$position <- as.integer(iam_leg2$position)

# See which positions in IAM leg are also in og leg
summary(iam_leg2$position %in% leg$position) # 470

# Remove rows where IAM position does not appear in leg
iam_leg3 <- iam_leg2[which(iam_leg2$position %in% leg$position), ]

# Remove rows of hap that no longer appear in IAM leg
iam_hap2 <- iam_hap[c(iam_leg3$row), ]

# Create new hap of zeros to add IAM info back into
iam_hap_out <- as.data.frame(matrix(0, nrow = nrow(leg), ncol = ncol(iam_hap2)))

# Merge IAM leg and og leg to see which rows to add IAM hap data to hap_out
leg_merge <- merge(iam_leg3, leg, by = 'position')

for (i in 1:nrow(leg_merge)) {
  
  # Add IAM data to rows where IAM leg positions overlap with leg positions
  iam_hap_out[leg_merge$row.y[i], ] <- iam_hap2[i, ]
}

write.table(hap1, paste0(dir, 'IAM_chr19_block17_coding_region.hap'), row.names = FALSE,
            col.names = FALSE, quote = FALSE, sep = '\t')


## Do the same thing for the IAM legend file

# Make empty legend file
iam_leg_out <- as.data.frame(matrix(nrow = nrow(leg), ncol = 4))
colnames(iam_leg_out) <- c('id', 'position', 'a0', 'a1')

# Add in necessray columns
iam_leg_out$position <- leg$position
iam_leg_out$a0 <- '.'
iam_leg_out$a1 <- '.'
iam_leg_out$id <- paste0('19:', iam_leg_out$position, '_Un_Known')

for(i in 1:nrow(leg_merge)){
  
  # Add IAM data to rows where IAM leg positions overlap with leg positions
  iam_leg_out[leg_merge$row.y[i],] <- leg_merge[i, c(2, 1, 3, 4)]
}

write.table(iam_leg_out, paste0(dir, 'IAM_chr19_block17_coding_region.legend'), 
            row.names = FALSE, quote = FALSE, sep = '\t')
