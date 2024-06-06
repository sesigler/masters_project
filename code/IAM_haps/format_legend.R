# This code is broken up into 2 parts
# Part 1: Format IAM legend file so it's compatible with hg19 block 37 leg files
# Step 1.1: Extract the IAM (hg38) coordinates and format them to be input into LiftOver 
#           to convert them to hg19 coordinates
# Step 1.2: Update the IAM legend file to have the hg19 coordinates
# Step 1.3: Remove positions from IAM leg and hap file not seen in the original 
#           legend file (using NFE_Block37_CDS_ref_added.legend from 
#           /storage/math/projects/compinfo/simulations/input/NFE_blocks since the
#           coding bp do not vary by ancestry).
#           Then format IAM hap and leg files so they include all the positions seen
#           in original legend file (add in filler data at the positions not seen in
#           the IAM data).
# Part 2:   Create IAM master legend file
# Step 2.1: Update unknown IAM variant positions with gnomAD AMR data
# Step 2.2: Update remaining unknown positions with annovar data
# Step 2.3: Add the functional annotation for each variant


# load libraries
library(dplyr)
library(tidyr)
library(data.table)

# main directory
dir = '/home/math/siglersa/IAM_haps/master_legend/'

# source of input IAM hap and leg files
dir_fin = '/home/math/siglersa/IAM_haps/final/'

# reference coding positions for chr19 block 37
dir_ref = '/storage/math/projects/compinfo/simulations/input/NFE_blocks/'

# AMR gnomAD data
dir_gnom = '/home/math/siglersa/IAM_haps/gnomad/'

# annovar functional annotations
dir_anno = '/home/math/siglersa/IAM_haps/annovar/'


# dir = dir_fin = dir_ref = dir_gnom = dir_anno = 'C:/Users/sagee/Documents/HendricksLab/IAM_hap_data/'


### Part 1: Subset IAM legend and hap files to all coding regions within block 37 (using code Megan Null gave me)


### Step 1.1: Write out positions of IAM legend file in LiftOver format for conversion to hg19 coordinates

# Read in IAM legend file
IAM_leg = read.table(paste0(dir_fin, 'hgdp_1kg_phased_haps_v2_block37_IAM.legend'), header=T, sep='')

# Create column of base pair (bp) positions in correct format for LiftOver
IAM_leg = IAM_leg %>% mutate(pos2 = paste0("chr19:", position, "-", position))

# extract formatted position column
out <- IAM_leg %>% select(pos2)

# save positions as txt file
write.table(out, paste0(dir, "IAM_legend_formatted_positons.txt"), quote=F, row.names=F, col.names=F)


### Step 1.2: Update IAM leg to have hg19 coordinates

# Read in IAM positions (hg38) converted to hg19 coordinates
iam_hg19 <- read.table(paste0(dir, 'IAM_leg_hg38_to_hg19_LiftOver_positions.bed'), header = F)

# Separate LiftOver positions into separate columns  
hg19_pos = iam_hg19 %>% separate(col="V1", into=c("chr", "pos1", "pos2"), sep = ":|-")

# Create new IAM leg with hg19 positions instead of hg38 coordinates
iam_leg_hg19 <- data.frame(position = hg19_pos$pos2, a0 = IAM_leg$a0, a1 = IAM_leg$a1)

# Make id column and make it first column
iam_leg_hg19 = iam_leg_hg19 %>% mutate(id = paste0("19:", position, "_", a0, "_", a1), 
                                       row = 1:nrow(iam_leg_hg19)) %>%
  select(id, everything())

# Save iam leg with hg19 coordinates
# write.table(iam_leg_hg19, paste0(dir, 'IAM_chr19_block37_hg19.legend'), 
#             row.names = FALSE, quote = FALSE, sep = '\t')


### Step 1.3: Add in missing coding bp and remove any bp not seen in original legend file 
###           to both IAM legend and hap files (Megan's code)

# Read in RAREsim legend file (chr19 block 37 positions are same across ancestry)
leg = read.table(paste0(dir_ref, 'NFE_Block37_CDS_ref_added.legend'), header=T, sep='\t')
leg$row <- 1:nrow(leg)

# Read in IAM leg file
# iam_leg = read.table(paste0(dir, 'IAM_chr19_block37_hg19.legend'), header=T, sep='\t')
iam_leg = iam_leg_hg19

# Read in IAM hap file
iam_hap = fread(paste0(dir_fin, 'hgdp_1kg_phased_haps_v2_block37_IAM.hap'))
iam_hap = as.data.frame(iam_hap)

# Check for duplicate positions in IAM legend file, then remove them
summary(duplicated(iam_leg$position)) # 131 (127 duplicates and 2 triplicates)
iam_leg2 <- iam_leg[-c(which(duplicated(iam_leg$position))), ]
iam_leg2$position <- as.integer(iam_leg2$position)

# See which positions in IAM leg are also in og leg
summary(iam_leg2$position %in% leg$position) # 470

# Remove rows where IAM position does not appear in leg
# Presumably these are non-coding positions
iam_leg3 <- iam_leg2[which(iam_leg2$position %in% leg$position), ]

# Remove rows of hap that no longer appear in IAM leg
iam_hap2 <- iam_hap[c(iam_leg3$row), ]

# Create new hap of zeros to add IAM info back into
iam_hap_out <- as.data.frame(matrix(0, nrow = nrow(leg), ncol = ncol(iam_hap2)))

# Merge IAM leg and og leg to see which rows to add IAM hap data to hap_out
leg_merge <- merge(iam_leg3, leg, by = 'position', suffixes = c("_iam", "_nfe"))

for (i in 1:nrow(leg_merge)) {
  
  # Add IAM data back into hap_out based on nfe row index
  iam_hap_out[leg_merge$row_nfe[i], ] <- iam_hap2[i, ]
}

# Save hap_out
write.table(iam_hap_out, paste0(dir, 'IAM_chr19_block37_coding_region.hap'), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = '\t')


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
  
  # Add IAM data back into leg_out based on nfe row index
  iam_leg_out[leg_merge$row_nfe[i],] <- leg_merge[i, c(2, 1, 3, 4)]
}

# Save legend file
write.table(iam_leg_out, paste0(dir, 'IAM_chr19_block37_coding_region.legend'), 
            row.names = FALSE, quote = FALSE, sep = '\t')


### Part 2: Create the IAM master legend file using Jess's code (process_data.Rmd)

# Step 2.1: gnomAD 

# read in reference legend file (19,029 bp)
leg.ref = read.table(paste0(dir, "IAM_chr19_block37_coding_region.legend"), header = TRUE) 
pos.hgdp = leg.ref %>% filter(!grepl("Un_Known", id)) %>% select(position) # 470 SNPs

# read in gnomad data (9,394 variants)
gnomad = read.table(paste0(dir_gnom, "gnomad.exomes.r2.1.1.sites.19.block37_AMR.vcf.INFO"), sep='\t', header=T) %>% rename(position = POS)
names(gnomad)[5:8] = sapply(strsplit(names(gnomad)[5:8], "_", fixed=T), head, 1)
pos.gnomad = gnomad %>% filter(position %in% leg.ref$position) %>% distinct(position) # 4,447 positions

# merge reference legend with gnomad data (19,525 lines)
combined = merge(leg.ref, gnomad, by="position", all.x=T) %>% arrange(position)

# subset the multiallelic SNVs
dups = combined %>% group_by(position) %>% filter(n()>1) # 473 alleles (496 duplicates)

# remove duplicates/triplicates from known HGDP multiallelic SNVs
dups.known = dups %>% filter(!grepl("Un_Known", id), a1==ALT) %>% mutate(prob = "1")

# extract the positions of unknown HGDP multiallelic SNVs
dups.unknown = dups %>% filter(grepl("Un_Known", id))
dup.pos = levels(as.factor(dups.unknown$position))

out = c()

# loop through the unknown HGDP multiallelic SNVs to choose one if possible from the gnomAD data
for (i in dup.pos){
  
  temp = dups.unknown %>% filter(position==i)
  
  # choose the allele with the maximum allele count in the population
  result = temp %>% filter(AC==max(AC))
  
  # if only one max, prob=1; if two maxes, prob=0.5
  result$prob = ifelse(nrow(result)==1, "1", "0.5")
  
  # if three maxes, prob=. (will use transition/transversion probabilities)
  result$prob = ifelse(nrow(result)==3, ".", result$prob)
  
  out = rbind(out, result)
}

# subset the biallelic SNVs
singles = combined %>% group_by(position) %>% filter(n()==1) 

# set the probability of gnomad and known HGDP SNVs to 1, otherwise .
singles$prob = ifelse(is.na(singles$REF), ".", "1")
singles$prob = ifelse(!grepl("Un_Known", singles$id), "1", singles$prob)

# combine biallelic SNVS back with the known/unknown multiallelic SNVs
combined2 = union(singles, union(dups.known, out)) %>% arrange(position)


### Step 2.2: ANNOVAR

# Read in ANNOVAR file
# Note this file is not ancestry specific, Jess just labeled it as NFE since that's
# the ancestry she started with
annovar = read.table(paste0(dir_anno, "annovar.chr19.block37.NFE.filtered.txt"), sep="\t")

# unknown variants not in gnomad or HGDP
unknown = combined2 %>% filter(grepl("Un_Known", id), is.na(REF)) %>% distinct(position)

# subset annovar to just the unknown alternate alleles
annovar.un = annovar %>% filter(V2 %in% unknown$position) 

# merge annovar variants with 1000G/gnomad variants
combined3 = merge(combined2, annovar.un, by.x="position", by.y="V2", all=T) %>% select(-V1, -V3)

### create master legend using the combined datasets
master = combined3
master$CHROM = "19"

# unknown positions
rows.un = which(grepl("Un_Known", master$id))

# loop through the unknown variants
for (i in rows.un){
  
  if (!is.na(master$REF[i])){
    
    # change a0/a1 from HGDP to the REF/ALT alleles from gnomad
    master[i, "a0"] = master[i, "REF"]
    master[i, "a1"] = master[i, "ALT"]
    
  } else if (!is.na(master$V4[i])){
    
    # change a0/a1 from HGDP to the REF/ALT alleles from annovar
    master[i, "a0"] = master[i, "V4"]
    master[i, "a1"] = master[i, "V5"]
  }
  
  # remove the "Un_Known" designation from the id
  id = substr(master[i, "id"], 1, nchar(master[i, "id"])-8)
  
  # rename the id based on the new REF/ALT alleles
  master[i, "id"] = paste0(id, master[i, "a0"], "_", master[i, "a1"])
}

master$AC = ifelse(is.na(master$AC), ".", master$AC)

# create file for annovar functional annotation
master2 = master %>% select(CHROM, START=position, END=position, REF=a0, ALT=a1, AC, prob)

write.table(master2, paste0(dir, 'master.chr19.block37.IAM.txt'),
row.names=F, col.names=F, quote=F, sep='\t')


### Need to run functional_annotation.sh in /home/math/siglersa/IAM_haps/annovar at this point


### Step 2.3: Functional Annotation

# read in functional annotation files
anno = read.table(paste0(dir_anno, "master.chr19.block37.IAM.txt.variant_function"), sep='\t') %>% 
  select(position2 = V4, InEx = V1, gene = V2) # all positions

anno.exo = read.table(paste0(dir_anno, "master.chr19.block37.IAM.txt.exonic_variant_function"), sep='\t') %>% 
  select(position2 = V5, fun = V2) # only exonic positions

# determine which positions are intronic
introns = which(anno$InEx=="intronic")

# add the functional annotation of the exons to the file with all positions
anno$fun = "."
anno[-introns, "fun"] = anno.exo$fun

# merge the functional annotations with the master legend file
leg.master = cbind(master, anno) %>% select(position, id, a0, a1, AC, prob, InEx, gene, fun)

write.table(leg.master, paste0(dir, 'chr19.block37.IAM.master.legend'),
row.names=F, col.names=F, quote=F, sep='\t')

summary(as.factor(leg.master$prob))


