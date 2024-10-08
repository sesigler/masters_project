### Used to plot power results ###

read_results <- function(dir, calc, scen, sub_scen, maf) {
  res <- read.csv(paste0(dir, calc, "_all_gene_", scen, "_", sub_scen, "_maf", maf, ".csv"), header=T)
  return(res)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpattern)
library(ggpubr)
library(binom)
library(data.table) # for fread

calc = 'Power'
Pops = c('IAM', 'NFE', 'EAS', 'AFR') # c('AFR', 'NFE')
admx_props = c(47, 44, 5, 4) # c(80, 20)
scen = 's2'
int_admx = '75% IAM, 19% NFE, 3% EAS, 3% AFR' # '100% AFR
ext_admx = '47% IAM, 44% NFE, 5% EAS, 4% AFR' # '80% AFR, 20% NFE' 
sub_scens = c('default', "160v100v85", "160v100v90", "160v100v95", "160v100v100")
# sub_scen = 'default'
comp = 'ccPrune'
adj = 'adjNeff'
pcase = 160
int_prune = 100
# ext_prune = 80
Ncase = 2000
Nic = 2000
Ncc = 10000
Nref = c(2000, 2000, 2000, 2000)
maf = 0.001 
genes_power = c("ADGRE5", "ADGRE3", "TECR") # genes used for cases (power)
str_Nrefs = 'Nref IAM: 2000, Nref NFE: 2000, Nref EAS: 2000, Nref AFR: 2000' # 'Nref AFR: 2000, Nref NFE: 2000'



# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/', tolower(calc), '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', tolower(calc), '/')
dir_out = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/power_plots/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/')

# file_in = paste0(scen, "_", sub_scen, "_maf", maf, ".csv")
# file_out = paste0(scen, "_", sub_scen, "_", adj, "_maf", maf, '.jpg')
file_out = paste0(scen, "_", comp, "_", adj, "_maf", maf, '.jpg')

# For single simulation scenario
res = read.csv(paste0(dir, calc, "_all_gene_", file_in), header=T)
res = pivot_longer(res, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)

# Read in results, for multiple scenario comparisons
res1 = read_results(dir, calc, scen, sub_scen=sub_scens[1], maf)
res2 = read_results(dir, calc, scen, sub_scen=sub_scens[2], maf)
res3 = read_results(dir, calc, scen, sub_scen=sub_scens[3], maf)
res4 = read_results(dir, calc, scen, sub_scen=sub_scens[4], maf)
res5 = read_results(dir, calc, scen, sub_scen=sub_scens[5], maf)
res6 = read_results(dir, calc, scen, sub_scen=sub_scens[6], maf)

# Format for ggplot
res1 = pivot_longer(res1, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = "default") 

res2 = pivot_longer(res2, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = sub_scens[2])

res3 = pivot_longer(res3, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = sub_scens[3])

res4 = pivot_longer(res4, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = sub_scens[4])

res5 = pivot_longer(res5, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = sub_scens[5])

res6 = pivot_longer(res6, prox_ext:burden_int, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, Nint = sub_scens[6])

res = rbind(res1, res2, res3, res4, res5)




results2 = res %>% mutate(MACs = rep(c("Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff", 
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Unadjusted", "Unadjusted"), times=60)) #12, 36, 48, 60, 72

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

# results2 = results2 %>% mutate(Data = rep(c("External", "External", "External",
#                                             "External", "External", "External",
#                                             "External", "External", "External",
#                                             "Internal + External", "Internal + External", "Internal + External",
#                                             "Internal + External", "Internal + External", "Internal + External",
#                                             "Internal", "Internal", "Internal"), 12))
# results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))

results2$MAF = factor(results2$MAF)
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Adjusted Ncc", "Adjusted Neff"), 
                       labels = c("Unadjusted", "Adjusted Ncc", "Adjusted"))
results2$Gene = factor(results2$Gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
                                               "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))

### version 1 plot labels
results2$casePrune = factor(results2$casePrune, levels=c("120v100v80", "140v100v80", "default"), 
                       labels = c("Cases (Power): 120% Pruned", "Cases (Power): 140% Pruned", "Cases (Power): 160% Pruned"))

results2$ccPrune = factor(results2$ccPrune, levels=c("default", "160v100v85", "160v100v90", "160v100v95", "160v100v100"), 
                            labels = c("Common Controls: 80% Pruned", "Common Controls: 85% Pruned", "Common Controls: 90% Pruned", "Common Controls: 95% Pruned", "Common Controls: 100% Pruned"))

results2$lessIC = factor(results2$lessIC, levels = c("default", "Ncase2K_Nic1K_Ncc1K", "Ncase2K_Nic1K_Ncc2K", "Ncase2K_Nic1K_Ncc5K", "Ncase2K_Nic1K_Ncc10K", "Ncase2K_Nic500_Ncc5K"),
                         labels = c("Cases: 2000, Internal Controls: 2000, Common Controls: 10000", "Cases: 2000, Internal Controls: 1000, Common Controls: 1000", "Cases: 2000, Internal Controls: 1000, Common Controls: 2000",
                                    "Cases: 2000, Internal Controls: 1000, Common Controls: 5000", "Cases: 2000, Internal Controls: 1000, Common Controls: 10000", "Cases: 2000, Internal Controls: 500, Common Controls: 5000"))

results2$lessIC2 = factor(results2$lessIC2, levels = c("default", "Ncase2K_Nic500_Ncc10K", "Ncase5K_Nic1K_Ncc10K", "Ncase5K_Nic500_Ncc10K"),
                          labels = c("Cases: 2000, Internal Controls: 2000, Common Controls: 10000", "Cases: 2000, Internal Controls: 500, Common Controls: 10000", 
                                     "Cases: 5000, Internal Controls: 1000, Common Controls: 10000", "Cases: 5000, Internal Controls: 500, Common Controls: 10000"))

results2$Nint = factor(results2$Nint, levels=c("Ncase500_Nic500", "Ncase3000_Nic3000", "Ncase1000_Nic1000", "Ncase4000_Nic4000", "default", "Ncase5000_Nic5000"),
                       labels = c("Cases and Internal Controls: 500", "Cases and Internal Controls: 3000", "Cases and Internal Controls: 1000",
                                  "Cases and Internal Controls: 4000", "Cases and Internal Controls: 2000", "Cases and Internal Controls: 5000"))

results2$Ncc = factor(results2$Ncc, levels=c("default", "Ncc20000", "Ncc50000"), 
                          labels = c("Common Controls: 10000", "Common Controls: 20000", "Common Controls: 50000"))

results2$Nref = factor(results2$Nref, levels=c("NrefAFR704_NrefNFE642", "default", "NrefAFR5000_NrefNFE5000", "NrefAFR10000_NrefNFE10000"), 
                          labels = c("Nref AFR: 704, Nref NFE: 642", "Nref AFR: 2000, Nref NFE: 2000", "Nref AFR: 5000, Nref NFE: 5000", "Nref AFR: 10000, Nref NFE: 10000"))
results2$Nref = factor(results2$Nref, levels=c("NrefIAM47_NrefNFE642_NrefEAS787_NrefAFR704", "default", "NrefIAM5K_NrefNFE5K_NrefEAS5K_NrefAFR5K", "NrefIAM10K_NrefNFE10K_NrefEAS10K_NrefAFR10K"), 
                       labels = c("Nref IAM: 47, Nref NFE: 642, Nref EAS: 787, Nref AFR: 704", "Nref IAM: 2K, Nref NFE: 2K, Nref EAS: 2K, Nref AFR: 2K", "Nref IAM: 5K, Nref NFE: 5K, Nref EAS: 5K, Nref AFR: 5K", "Nref IAM: 10K, Nref NFE: 10K, Nref EAS: 10K, Nref AFR: 10K"))

### version 2 plot labels
results2$casePrune = factor(results2$casePrune, levels=c("120v100v80", "140v100v80", "default"), 
                            labels = c("120%", "140%", "160%"))

results2$ccPrune = factor(results2$ccPrune, levels=c("default", "160v100v85", "160v100v90", "160v100v95", "160v100v100"), 
                          labels = c("80%", "85%", "90%", "95%", "100%"))

results2$lessIC = factor(results2$lessIC, levels = c("default", "Ncase2K_Nic1K_Ncc1K", "Ncase2K_Nic1K_Ncc2K", "Ncase2K_Nic1K_Ncc5K", "Ncase2K_Nic1K_Ncc10K", "Ncase2K_Nic500_Ncc5K"),
                         labels = c("(2K, 2K, 10K)", "(2K, 1K, 1K)", "(2K, 1K, 2K)", "(2K, 1K, 5K)", "(2K, 1K, 10K)", "(2K, 500, 5K)"))

results2$lessIC2 = factor(results2$lessIC2, levels = c("default", "Ncase2K_Nic500_Ncc10K", "Ncase5K_Nic1K_Ncc10K", "Ncase5K_Nic500_Ncc10K"),
                          labels = c("(2K, 2K, 10K)", "(2K, 500, 10K)", "(5K, 1K, 10K)", "(5K, 500, 10K)"))

results2$Nint = factor(results2$Nint, levels=c("Ncase500_Nic500", "Ncase1000_Nic1000", "default", "Ncase3000_Nic3000", "Ncase4000_Nic4000", "Ncase5000_Nic5000"),
                       labels = c("500", "1K", "2K", "3K", "4K", "5K"))

results2$Ncc = factor(results2$Ncc, levels=c("default", "Ncc20000", "Ncc50000"), 
                      labels = c("10K", "20K", "50K"))

results2$Nref = factor(results2$Nref, levels=c("NrefAFR704_NrefNFE642", "default", "NrefAFR5000_NrefNFE5000", "NrefAFR10000_NrefNFE10000"), 
                       labels = c("(704, 642)", "(2000, 2000)", "(5000, 5000)", "(10000, 10000)"))
results2$Nref = factor(results2$Nref, levels=c("NrefIAM47_NrefNFE642_NrefEAS787_NrefAFR704", "default", "NrefIAM5K_NrefNFE5K_NrefEAS5K_NrefAFR5K", "NrefIAM10K_NrefNFE10K_NrefEAS10K_NrefAFR10K"), 
                       labels = c("(47, 642, 787, 704)", "(2K, 2K, 2K, 2K)", "(5K, 5K, 5K, 5K)", "(10K, 10K, 10K, 10K)"))


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
pdefault <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
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
  facet_wrap(~Data, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene \nCases: ', pcase,  '% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune,  '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref AFR: ', Nref[1], ', Nref NFE: ', Nref[2], ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pdefault
ggsave(file = paste0(dir_out, calc, '_gene_', file_out), plot = pdefault, height = 8, width = 16, units = 'in')

### Version 1 plots (Genes on x-axis, facet_wrap by comp)

# power for new admixed AFR results-Cases (Power) % Pruned
pcasePrune <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                     aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~casePrune, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Cases (Power) Pruning \nCases: ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                                '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pcasePrune
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = pcasePrune, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Common Controls % Pruned
pccPrune <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                     aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~ccPrune, ncol = 2, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Common Control Pruning \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pccPrune
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = pccPrune, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Less Internal Controls than Common Controls
plessIC <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~lessIC, ncol = 2, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Case, Internal Control, and Common Control Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\n', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
plessIC
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = plessIC, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Less Internal Controls than Common Controls v2
plessIC2 <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                  aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~lessIC2, ncol = 2, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Case, Internal Control, and Common Control Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\n', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
plessIC2
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = plessIC2, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Internal Sample Size
pNint <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                   aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Nint, ncol = 2, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Internal Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcc: ', Ncc, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pNint
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = pNint, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Common Control Sample Size
pNcc <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                   aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Ncc, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Common Control Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pNcc
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = pNcc, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Reference Sample Size
pNref <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
               aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Nref, ncol = 2, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Gene', title=paste0('Power by Gene and Reference Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pNref
ggsave(file = paste0(dir_out, calc, '_v1_gene_', file_out), plot = pNref, height = 8, width = 16, units = 'in')


### Version 2 plots (Comp on x-axis, facet_wrap by gene)

# power for new admixed AFR results-Cases (Power) % Pruned
pcasePrune <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                     aes(x=casePrune, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Cases (Power): % Pruned', title=paste0('Power by Gene and Cases (Power) Pruning \nCases: ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pcasePrune
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = pcasePrune, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Common Controls % Pruned
pccPrune <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                   aes(x=ccPrune, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Common Controls: % Pruned', title=paste0('Power by Gene and Common Control Pruning \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pccPrune
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = pccPrune, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Less Internal Controls than Common Controls
plessIC <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
               aes(x=lessIC, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Sample Size (Cases, Internal Controls, Common Controls)', title=paste0('Power by Gene and Case, Internal Control, and Common Control Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                                               '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\n', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
plessIC
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = plessIC, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Less Internal Controls than Common Controls v2
plessIC2 <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                  aes(x=lessIC2, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Sample Size (Cases, Internal Controls, Common Controls)', title=paste0('Power by Gene and Case, Internal Control, and Common Control Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                                                                            '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\n', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
plessIC2
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = plessIC2, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Internal Sample Size
pNint <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                aes(x=Nint, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Internal Sample Size', title=paste0('Power by Gene and Internal Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcc: ', Ncc, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pNint
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = pNint, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Common Control Sample Size
pNcc <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
               aes(x=Ncc, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Common Control Sample Size', title=paste0('Power by Gene and Common Control Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', ', str_Nrefs, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pNcc
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = pNcc, height = 8, width = 16, units = 'in')

# power for new admixed AFR results-Reference Sample Size
pNref <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc") & (Gene == "ADGRE5" | Gene == "ADGRE3" | Gene == "TECR")), 
                aes(x=Nref, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.80, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~Gene, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Power', x='Reference Sample Size (AFR, NFE)', title=paste0('Power by Gene and Reference Sample Size \nCases: 160% pruned, ', int_admx, '\nInternal Controls: ', int_prune, '% pruned, ', int_admx,
                                         '\nCommon Controls: ', ext_prune, '% pruned, ', ext_admx, '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
pNref
ggsave(file = paste0(dir_out, calc, '_v2_gene_', file_out), plot = pNref, height = 8, width = 16, units = 'in')

