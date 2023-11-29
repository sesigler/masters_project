# This file is used to extract information about each gene in
# chromosome 19 Block 37, such as number of variants (common and rare),
# number of fun/syn variants, etc.

# load libraries
library(dplyr)
library(tidyr)
library(data.table)

Pop = 'NFE'
p_exp = 100
i = 1
maf = 0.01

dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/gene_info/')

leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.', p_exp, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t')
leg$row = 1:nrow(leg)

# meg_leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.copy.legend'), header=T, sep='\t')
# leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene) %>% count()
# meg_leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene) %>% count()

hap = fread(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.all.', p_exp, 'fun.', p_exp, 'syn.haps.gz'))
hap = as.data.frame(hap)

# add relevant columns to a copy of the hap file
hap2 = hap %>% mutate(sums = rowSums(hap), af = sums/ncol(hap), gene = leg$gene, fun = leg$fun, row = leg$row)

# See how many variants there are by gene
hap2 %>% filter(sums != 0) %>% group_by(gene) %>% count()

# See how many common variants there are by gene
hap2 %>% filter(af > maf) %>% group_by(gene) %>% count()

# See how many fun/syn variants there are by gene
fun = hap2 %>% group_by(gene, fun) %>% count()
