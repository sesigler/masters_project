# This code looks into the number of variants used for iECAT and the corresponding
# p-value

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

Pop1 = "AFR"
Pop2 = "NFE"
# scen = "s1"
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc10k'  #Number of common controls: 'cc5k' or 'cc10k'
int_prune = 100
ext_prune = 80
pruning = "pruneSepRaresim" #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '120v100v80'
data = 'by_gene'

dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/')
iecat_info = read.table(paste0(dir, "iecat_info_", folder, "_", Pop2, "_maf", maf, ".txt"), header = T)

plot(iecat_info$n.marker.test, iecat_info$p.value)

# Scatter plot of number of variants used v p-value by gene
ggplot(iecat_info, 
       aes(x=n.marker.test, y=p.value, color=Gene)) +
  geom_point() +
  facet_wrap(~Gene) +
  labs(y='P-value', x='Number of Variants Used', title=paste0('Number of Variants Used and Corresponding P-value for iECAT')) +
  theme_bw(base_size = 15)

# Bar plot of number of variants used by gene
ggplot(iecat_info, 
       aes(x=n.marker.test, color=Gene)) +
  geom_bar(stat = 'bin', fill='white') +
  facet_wrap(~Gene) +
  labs(y='Count', x='Number of Variants Used', title=paste0('Number of Variants Used and Corresponding P-value for iECAT')) +
  theme_bw(base_size = 15)

# Density plot of number of variants used by gene
ggplot(iecat_info, 
       aes(x=n.marker.test, color=Gene)) +
  geom_density() +
  facet_wrap(~Gene) +
  labs(y='Density', x='Number of Variants Used', title=paste0('Density of Number of Variants Used by Gene for iECAT')) +
  theme_bw(base_size = 15)
