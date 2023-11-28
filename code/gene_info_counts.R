# load libraries
library(dplyr)
library(tidyr)
library(data.table)

Pop = 'AFR'
p_exp = 100
i = 1

dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/gene_info/')

leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.', p_exp, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t')
leg$row = 1:nrow(leg)

meg_leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.copy.legend'), header=T, sep='\t')

hap = fread(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.all.', p_exp, 'fun.', p_exp, 'syn.haps.gz'))
hap = as.data.frame(hap)

hap2 = hap
hap2$sums = rowSums(hap2)
hap2$row = leg$row
hap2$gene = leg$gene

hap2 %>% filter(sums != 0) %>% group_by(gene) %>% count()

leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene) %>% count()
meg_leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene) %>% count()
