# This file is used to extract information about each gene in
# chromosome 19 Block 37, such as number of variants (common and rare),
# number of fun/syn variants, etc.

# load libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)

Pop = 'AFR'
Pop1 = 'AFR'
Pop2 = 'NFE'
p_exp = 100
i = 1
maf = 0.01

dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/gene_info/')

### RAREsim variants
leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.', p_exp, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t')
leg$row = 1:nrow(leg)

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

# Check which RAREsim variants are unique to AFR and NFE pops
leg_afr = read.table(paste0(dir, 'chr19.block37.', Pop1, '.sim', i, '.', p_exp, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t')
leg_afr$row = 1:nrow(leg_afr)

leg_nfe = read.table(paste0(dir, 'chr19.block37.', Pop2, '.sim', i, '.', p_exp, 'fun.', p_exp, 'syn.legend'), header=T, sep='\t')
leg_nfe$row = 1:nrow(leg_nfe)

hap_afr = fread(paste0(dir, 'chr19.block37.', Pop1, '.sim', i, '.all.', p_exp, 'fun.', p_exp, 'syn.haps.gz'))
hap_afr = as.data.frame(hap_afr)

hap_nfe = fread(paste0(dir, 'chr19.block37.', Pop2, '.sim', i, '.all.', p_exp, 'fun.', p_exp, 'syn.haps.gz'))
hap_nfe = as.data.frame(hap_nfe)

hap_afr = hap_afr %>% mutate(sums = rowSums(hap_afr), af = sums/ncol(hap_afr), gene = leg_afr$gene, fun = leg_afr$fun, row = leg_afr$row, position = leg_afr$position)
hap_nfe = hap_nfe %>% mutate(sums = rowSums(hap_nfe), af = sums/ncol(hap_nfe), gene = leg_nfe$gene, fun = leg_nfe$fun, row = leg_nfe$row, position = leg_nfe$position)

afr_vars = hap_afr %>% filter(sums != 0)
nfe_vars = hap_nfe %>% filter(sums != 0)

afr_uniq = setdiff(afr_vars$position, nfe_vars$position)
nfe_uniq = setdiff(nfe_vars$position, afr_vars$position)

afr_variants = afr_vars %>% filter(position %in% afr_uniq)
nfe_variants = nfe_vars %>% filter(position %in% nfe_uniq)

tst1 = afr_variants %>% group_by(gene) %>% count()
tst2 = nfe_variants %>% group_by(gene) %>% count()

tst1 = tst1 %>% mutate(Pop = 'AFR')
tst2 = tst2 %>% mutate(Pop = 'NFE')
tst = rbind(tst1, tst2)

ggplot(data=tst, aes(x=gene, y=n, fill=Pop)) +
  geom_bar(stat="identity", position=position_dodge())

### gnomAD variants
meg_leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', i, '.copy.legend'), header=T, sep='\t')
meg_leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene) %>% count()
meg_leg_fun = meg_leg %>% filter(AC!='.' & AC!=0) %>% group_by(gene, fun) %>% count()

# Check which gromAD variants are unique to AFR and NFE pops
meg_leg_afr = read.table(paste0(dir, 'chr19.block37.', Pop1, '.sim', i, '.copy.legend'), header=T, sep='\t')
meg_leg_afr$row = 1:nrow(meg_leg_afr)

meg_leg_nfe = read.table(paste0(dir, 'chr19.block37.', Pop2, '.sim', i, '.copy.legend'), header=T, sep='\t')
meg_leg_nfe$row = 1:nrow(meg_leg_nfe)

ml_afr = meg_leg_afr %>% filter(AC!='.' & AC!=0)
ml_nfe = meg_leg_nfe %>% filter(AC!='.' & AC!=0)

afr_vars = setdiff(ml_afr$row, ml_nfe$row)
nfe_vars = setdiff(ml_nfe$row, ml_afr$row)

afr_variants = ml_afr %>% filter(row %in% afr_vars)
nfe_variants = ml_nfe %>% filter(row %in% nfe_vars)

tst1 = afr_variants %>% group_by(gene) %>% count()
tst2 = nfe_variants %>% group_by(gene) %>% count()

tst1 = tst1 %>% mutate(Pop = 'AFR')
tst2 = tst2 %>% mutate(Pop = 'NFE')
tst = rbind(tst1, tst2)

ggplot(data=tst, aes(x=gene, y=n, fill=Pop)) +
  geom_bar(stat="identity", position=position_dodge())
