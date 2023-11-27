# load libraries
library(dplyr)
library(tidyr)
library(data.table)

Pop = 'AFR'
p_exp = 100
i = 1

dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/gene_info/')

leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.', p_exp, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t')

hap = fread(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.all.', p_exp, 'fun.', p_exp, 'syn.haps.gz'))
hap = as.data.frame(hap)

leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene) %>% count()
