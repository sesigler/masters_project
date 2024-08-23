# For plotting the relationships between the fun:syn ratios in cases and controls and the p-value

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(binom)
library(data.table)

Pop_admx = 'AFR_NFE'  
Pops = c('AFR', 'NFE')
admx_props = c(80, 20)
scen = 's2'
sub_scen = 'default'
maf = 0.001 

dir = paste0('C:/Users/sagee/Documents/HendricksLab/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/')

acs = read.csv(paste0(dir, "ACs_prox_ext_adj_Ncc_", scen, "_", sub_scen, "_maf", maf, ".csv"), header=T)

acs2 <- acs %>% mutate(case_control_ratio = (case_fun/case_syn)/(control_fun/control_syn))

p1 <- ggplot(acs2 %>% filter((gene == "ADGRE5" | gene == "ADGRE3" | gene == "TECR")), 
             aes(x=prox, y=case_control_ratio)) +
  geom_point() +
  geom_vline(xintercept=0.05, linetype=2, linewidth=0.5) +
  facet_wrap(~gene, ncol = 1) +
  labs(y='Case:Control Fun Syn Ratio', x='P-value') +
  theme_bw(base_size = 14)
p1

p2 <- ggplot(acs2 %>% filter((gene == "ADGRE5" | gene == "ADGRE3" | gene == "TECR")), 
             aes(x=control_ratio, y=case_ratio)) +
  geom_point() +
  geom_abline() +
  facet_wrap(~gene, ncol = 3) +
  labs(y='Case Ratio', x='Control Ratio') +
  theme_bw(base_size = 14)
p2
