# Create a list of the loci that are non-zero in the reference haplotypes to 
# use as the disease loci for HAPGEN2

# load libraries
library(dplyr)
library(tidyr)
library(data.table)

dir = '/home/math/siglersa/IAM_haps/master_legend/'
# dir = 'C:/Users/sagee/Documents/HendricksLab/IAM_hap_data/'

# Read in reference hap file
hap = fread(paste0(dir, 'IAM_chr19_block37_coding_region.hap'))

# Read in legend file
leg = read.table(paste0(dir, 'IAM_chr19_block37_coding_region.legend'), header=T, sep='\t')

# Sum up all the alternate alleles in each row
loci = data.frame(position = leg$position, sums=rowSums(hap))

# Select only the non-zero rows
loci2 = loci %>% filter(sums != 0)

# Write out non-zero loci
write.table(loci2, paste0(dir, "IAM_disease_loci.txt"), quote=F, row.names=F, col.names=T)