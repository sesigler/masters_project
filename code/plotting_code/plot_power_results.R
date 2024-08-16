### Used to plot power results ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(binom)
library(data.table) # for fread

calc = 'Power'
Pops = c('AFR', 'NFE')
admx_props = c(80, 20)
scen = 's1'
sub_scens = c('default', '160v100v85', '160v100v90', '160v100v95')
sub_scen = 'default'
comp = 'ccPrune'
adj = 'adjNeff'
pcase = 160
int_prune = 100
ext_prune = 80
Ncase = Nic = 2000
Ncc = 10000
Nref = c(2000, 2000)
maf = 0.001 
genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)



# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/', tolower(calc), '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', tolower(calc), '/')
dir_out = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/power_plots/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/')

file_in = paste0(scen, "_", sub_scen, "_maf", maf, ".csv")
file_out = paste0(scen, "_", sub_scen, "_", adj, "_maf", maf, '.jpg')
file_out = paste0(scen, "_", comp, "_", adj, "_maf", maf, '.jpg')


res = read.csv(paste0(dir, calc, "_all_gene_", file_in), header=T)
res = pivot_longer(res, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)



t1e_gene1 = pivot_longer(t1e_gene1, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene2 = pivot_longer(t1e_gene2, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene3 = pivot_longer(t1e_gene3, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
# t1e_gene = pivot_longer(t1e_gene, prox_int:burden_all, names_to="Method", values_to="Value") %>% 
#   mutate(MAF = maf)

t1e_gene = rbind(t1e_gene1, t1e_gene2, t1e_gene3)

results2 = res %>% mutate(MACs = rep(c("Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff", 
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Unadjusted", "Unadjusted"), times=12))

results2$Method = factor(results2$Method, levels=c("prox_ext", "prox_ext_adj_Ncc", "prox_ext_adj_Neff", 
                                                   "proxW_ext", "proxW_ext_adj_Ncc", "proxW_ext_adj_Neff",
                                                   "prox2_ext", "prox2_ext_adj_Ncc", "prox2_ext_adj_Neff", 
                                                   "prox2_all", "prox2_all_adj_Ncc", "prox2_all_adj_Neff",
                                                   "iecat_all", "iecat_all_adj_Ncc", "iecat_all_adj_Neff",
                                                   "skato_int", "skat_int", "burden_int"),
                         labels=rep(c("ProxECAT (External)", "ProxECAT-weighted (External)", 
                                      "LogProx (External)", "LogProx (Internal + External)", 
                                      "iECAT-O (Internal + External)", 
                                      "SKAT-O (Internal)", "SKAT (Internal)", "Burden (Internal)"), 
                                    times=c(3, 3, 3, 3, 3, 1, 1, 1)))

# results2 = results2 %>% mutate(Data = rep(c("Internal", "External", "External", "External",
#                                             "Internal", "External", "External", "External",
#                                             "Internal", "External", "Internal + External", "External", "External", "Internal + External", "Internal + External",
#                                             "Internal + External", "Internal + External", "Internal + External",
#                                             "Internal", "External", "Internal + External",
#                                             "Internal", "External", "Internal + External",
#                                             "Internal", "External", "Internal + External"), 36))
# results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))


results2 = results2 %>% mutate(Nref = rep(c("AFR: 704, NFE: 642", "AFR: 2000, NFE: 2000", "AFR: 10000, NFE: 10000"), each=324))

results2$MAF = factor(results2$MAF)
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Adjusted Ncc", "Adjusted Neff"), 
                       labels = c("Unadjusted", "Adjusted Ncc", "Adjusted"))
results2$Gene = factor(results2$Gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
                                               "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))


# get CI's
results2$Lower = '.'
results2$Upper = '.'

## CI stack exchange:
# https://stats.stackexchange.com/questions/82720/confidence-interval-around-binomial-estimate-of-0-or-1
# paper: https://projecteuclid.org/journals/statistical-science/volume-16/issue-2/Interval-Estimation-for-a-Binomial-Proportion/10.1214/ss/1009213286.full?tab=ArticleFirstPage

# default level is 95% confidence
nsim = 100
for(i in 1:nrow(results2)){
  results2$Lower[i] = binom.confint(nsim*results2$Value[i], nsim, method=c("wilson"), type="central")$lower
  results2$Upper[i] = binom.confint(nsim*results2$Value[i], nsim, method=c("wilson"), type="central")$upper
}

results2$Lower = as.numeric(results2$Lower)
results2$Upper = as.numeric(results2$Upper)

# colorblind friendly palette
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors2 = c("#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#CC79A7", "#999999", "#BC9F4C", "#009E73")

# Controls everything in the graph
# base_size = 25


# t1e for new admixed AFR results
p1 <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
              aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.8, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  # facet_wrap(~ccPrune, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene (Default Parameters) \nCases: ', pcase,  '% pruned, 80% AFR, 20% NFE\nInternal Controls: ', int_prune, '% pruned, 80% AFR, 20% NFE',
                                                '\nCommon Controls: ', ext_prune,  '% pruned, 80% AFR, 20% NFE\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref AFR: ', Nref[1], ', Nref NFE: ', Nref[2], ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
p1
ggsave(file = paste0(dir_out, calc, '_gene_', file_out), plot = p1, height = 8, width = 16, units = 'in')
