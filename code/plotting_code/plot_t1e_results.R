### Used to plot type I error results ###

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

calc = 'T1e'
Pops = c('AFR', 'NFE')
admx_props = c(80, 20)
scen = 's2'
sub_scens = c('default', '160v100v85', '160v100v90', '160v100v95')
comp = 'ccPrune'
adj = 'adjNeff'
int_prune = 100
# ext_prune = 80
Ncase = Nic = 2000
Ncc = 10000
Nref = c(2000, 2000)
maf = 0.001 


# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/', tolower(calc), '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', tolower(calc), '/')
dir_out = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/typeI_error_plots/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/')

# file_in = paste0(scen, "_", sub_scen, "_maf", maf, ".csv")
file_out = paste0(scen, "_", comp, "_", adj, "_maf", maf, '.jpg')

# read in the results
# Proportion Estimates 
# t1e_case_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_case_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_cc_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_cc_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_int_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_int_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

# results = read.csv(paste0(dir, calc, "_all_gene_", file_in), header=T)
# results = pivot_longer(results, prox_int:burden_all, names_to="Method", values_to="Value") %>%
#   mutate(MAF = maf)

# Read in results
res1 = read_results(dir, calc, scen, sub_scen=sub_scens[1], maf)
res2 = read_results(dir, calc, scen, sub_scen=sub_scens[2], maf)
res3 = read_results(dir, calc, scen, sub_scen=sub_scens[3], maf)
res4 = read_results(dir, calc, scen, sub_scen=sub_scens[4], maf)

# Format for ggplot
res1 = pivot_longer(res1, prox_ext:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = "160v100v80") 

res2 = pivot_longer(res2, prox_ext:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = "160v100v85")

res3 = pivot_longer(res3, prox_ext:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = "160v100v90")

res4 = pivot_longer(res4, prox_ext:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf, ccPrune = "160v100v95")

res = rbind(res1, res2, res3, res4)


results2 = res %>% mutate(MACs = rep(c("Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff", 
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Adjusted Ncc", "Adjusted Neff",
                                       "Unadjusted", "Unadjusted", "Unadjusted",
                                       "Unadjusted", "Unadjusted", "Unadjusted",
                                       "Unadjusted", "Unadjusted", "Unadjusted"), times=48)) #36 or 48

results2$Method = factor(results2$Method, levels=c("prox_ext", "prox_ext_adj_Ncc", "prox_ext_adj_Neff", 
                                                   "proxW_ext", "proxW_ext_adj_Ncc", "proxW_ext_adj_Neff",
                                                   "prox2_ext", "prox2_ext_adj_Ncc", "prox2_ext_adj_Neff", 
                                                   "prox2_all", "prox2_all_adj_Ncc", "prox2_all_adj_Neff",
                                                   "iecat_all", "iecat_all_adj_Ncc", "iecat_all_adj_Neff",
                                                   "skato_int", "skato_ext", "skato_all",
                                                   "skat_int", "skat_ext", "skat_all",
                                                   "burden_int", "burden_ext", "burden_all"),
                         labels=rep(c("ProxECAT (External)", "ProxECAT-weighted (External)", 
                                      "LogProx (External)", "LogProx (Internal + External)", 
                                      "iECAT-O (Internal + External)", 
                                      "SKAT-O (Internal)", "SKAT-O (External)", "SKAT-O (Internal + External)",
                                      "SKAT (Internal)", "SKAT (External)", "SKAT (Internal + External)", 
                                      "Burden (Internal)", "Burden (External)", "Burden (Internal + External)"), times=c(3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1)))

# results2 = results2 %>% mutate(Data = rep(c("External", "External", "External",
#                                             "External", "External", "External",
#                                             "External", "External", "External", 
#                                             "Internal + External", "Internal + External", "Internal + External",
#                                             "Internal + External", "Internal + External", "Internal + External",
#                                             "Internal", "External", "Internal + External",
#                                             "Internal", "External", "Internal + External",
#                                             "Internal", "External", "Internal + External"), 36))
# results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))

results2$MAF = factor(results2$MAF)
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Adjusted Ncc", "Adjusted Neff"), 
                       labels = c("Unadjusted", "Adjusted Ncc", "Adjusted"))
results2$Gene = factor(results2$Gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
                                               "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))
results2$ccPrune = factor(results2$ccPrune, levels=c("160v100v80", "160v100v85", "160v100v90", '160v100v95'), 
                          labels = c("Common Controls: 80% Pruned", "Common Controls: 85% Pruned", "Common Controls: 90% Pruned", "Common Controls: 95% Pruned"))

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
cbPalette = c("#999999", "#BC9F4C", "#CC79A7", "#D55E00", "#E69F00", "#F0E442", "#009E73", "#56B4E9", "#0072B2")

colors_nfe = c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999",
               "lightgreen", "#56B4E9", "#E69F00", "pink")

colors_meth = c("#56B4E9", "#0072B2", "lightgreen", "#009E73", "#CC79A7", "#BC9F4C", "#E69F00", "#D55E00")

colors2 = c("#56B4E9", "#0072B2", "#E69F00", "#D55E00", "#CC79A7", "#999999", "#BC9F4C", "#009E73")

colors3 = c("#56B4E9", "#009E73", "#E69F00", "#D55E00", "#CC79A7", "#999999", "#BC9F4C", "#0072B2")


colors_gene = c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999",
                "lightgreen", "#56B4E9", "#E69F00", "#BC9F4C", "pink", "red", "purple")

 

# t1e for admixed by gene results from pcasev100vpconf pipeline-a for all methods
p10a <- ggplot(results2 %>% filter(MACs != paste0("Variant Adjusted Ncc ", Ncc) & ((Method == "ProxECAT" & Data == "External") | (Method == "ProxECAT-weighted" & Data == "External") | 
                                     (Method == "LogProx" & (Data == "External" | Data == "Internal + External")) |
                                     (Method == "iECAT-O") |
                                     ((Method == "SKAT-O" | Method == "SKAT" | Method == "Burden") & Data == "Internal"))), 
               aes(x=Gene, y=Value, color=Method, shape=MACs)) +
          geom_point(size=2, position=position_dodge(width=0.5)) +
          geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
          # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
          # scale_y_continuous(limits=c(0, 1)) +
          # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
          # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
          geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
          scale_color_manual(values=colors_meth) +
          scale_shape_manual(values = c(1, 16)) +
          # scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
          # scale_size_manual(values = c(0.5, 1, 3)) +
          facet_wrap(~Data, ncol = 1, scales = 'free_y') +
          # facet_wrap(~Data, ncol = 1) +
          labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Controls Used \nCases and Internal Controls: ', int_prune, '% pruned, 100% AFR \nCommon Controls: ', 
                                                        ext_prune, '% pruned, 80% AFR, 20% NFE', 
                                                        '\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: ', maf)) +
          # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
          theme_bw(base_size = 15)
p10a
ggsave(file = paste0(dir_out, calc, '_', file_out), plot = p10a, height = 8, width = 16, units = 'in')

# t1e for admixed by gene results from pcasev100vpconf pipeline-m for main methods
p10m <- ggplot(results2 %>% filter(!(Method == "SKAT-O" | Method =="SKAT" | Method == "Burden")), aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=colors_meth) +
  scale_shape_manual(values = c(9, 1, 16)) +
  # scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  # scale_size_manual(values = c(0.5, 1, 3)) +
  facet_wrap(~Data, ncol = 1, scales = 'free_y') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Controls Used \nInternal Data ', int_prune, '% pruned vs External Data ', ext_prune, '% pruned from ', folder, 
                                                ' Data \nScenario: ', scen, ', Pruning: ', pruning_plot, 
                                                '\nCases and Internal Controls: 80% AFR 20% NFE, External Controls: 80% AFR 20% NFE \nNsim: ', Nsim, ', Ncase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref ', Pop1, ': ', Nref_pop1, ', Nref ', Pop2, ': ', Nref_pop2, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 15)
p10m
ggsave(file = paste0(dir_out_t1e, file_out), plot = p10m, height = 8, width = 16, units = 'in')

# t1e for different adjustments from pcasev100vpconf pipeline-proxECAT only
p11 <- ggplot(results2 %>% filter(!(MACs == "Gene Adjusted Ncc" | MACs == "Gene Adjusted Neff")), aes(x=Gene, y=Value, color=MACs, shape=Method)) +
          geom_point(size=3, position=position_dodge(width=0.5)) +
          geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
          # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
          # scale_y_continuous(limits=c(0, 1)) +
          # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
          # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
          geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
          scale_color_manual(values=colors_prox) +
          scale_shape_manual(values = c(1, 16)) +
          facet_wrap(~Nref+Data, ncol = 2, scales = 'fixed', labeller = label_wrap_gen(multi_line=FALSE)) +
          # facet_wrap(~Data, ncol = 1) +
          labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene: ', int_prune, '% vs ', ext_prune, '% from ', folder, 
                                                        ' Data \nScenario: ', scen, ', Pruning: ', pruning_plot, 
                                                        '\nCases and Internal Controls: 80% AFR 20% NFE, External Controls: 80% AFR 20% NFE \nNsim: ', Nsim, ', Ncase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref: Varies, MAF: ', maf)) +
          theme(axis.text.x = element_text(angle = 35, hjust=0.65))
          # theme_bw(base_size = 15)
p11
ggsave(file = paste0(dir_out_t1e, file_out), plot = p11, height = 8, width = 16, units = 'in')

# t1e for admixed by gene results for STATGEN-EXTERNAL ONLY
p12 <- ggplot(results2 %>% filter(!(Method == "SKAT-O" | Method =="SKAT" | Method == "Burden")) %>% 
                filter(Data == "External") %>%
                filter(!(MACs == "Variant Adjusted Neff")), 
              aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=colors_meth) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  # scale_size_manual(values = c(0.5, 1, 3)) +
  facet_wrap(~Nref, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Reference Sample Size: Cases vs Common Controls \nCases: ', int_prune, '% pruned, 100% AFR \nCommon Controls: ', 
                                                ext_prune, '% pruned, 80% AFR, 20% NFE', 
                                                '\nNsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 15)
p12
ggsave(file = paste0(dir_out_t1e, 't1e_gene_Nref_ext_only_', file_out), plot = p12, height = 8, width = 16, units = 'in')

# t1e for admixed by gene results for STATGEN-INTERNAL+EXTERNAL ONLY
p13 <- ggplot(results2 %>% filter(!(Method == "SKAT-O" | Method =="SKAT" | Method == "Burden")) %>% 
                filter(Data == "Internal + External") %>%
                filter(!(MACs == "Variant Adjusted Neff")), 
              aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=3, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
  scale_color_manual(values=colors_meth) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  # scale_size_manual(values = c(0.5, 1, 3)) +
  facet_wrap(~Nint, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Internal Sample Size \nCases and Internal Controls: ', int_prune, '% pruned, 100% AFR \nCommon Controls: ', 
                                                ext_prune, '% pruned, 80% AFR, 20% NFE', 
                                                '\nNcc: ', Ncc, ', MAF: ', maf, ', Nref AFR: ', Nref, ', Nref NFE: ', Nref)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 15)
p13
ggsave(file = paste0(dir_out_t1e, 't1e_gene_Nref_int_and_ext_only_', file_out), plot = p13, height = 8, width = 16, units = 'in')

# t1e for new admixed AFR results
p14 <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc" | Method == "SKAT-O (External)" | Method == "SKAT-O (Internal + External)" |
                                      Method == "SKAT (External)" | Method == "SKAT (Internal + External)" |
                                      Method == "Burden (External)" | Method == "Burden (Internal + External)")), 
              aes(x=Gene, y=Value, color=Method, shape=MACs)) +
  geom_point(size=1.8, position=position_dodge(width=0.8)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=0.5) +
  # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
  # scale_y_continuous(limits=c(0, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
  # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=0.5, width=.5, position=position_dodge(width=0.8)) +
  scale_color_manual(values = colors2) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("dashed", "solid")) +
  # scale_size_manual(values = c(0.5, 0.8)) +
  facet_wrap(~ccPrune, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Common Control Pruning \nCases and Internal Controls: ', int_prune, '% pruned, 100% AFR',
                                                '\nCommon Controls: 80% AFR, 20% NFE\nNcase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref AFR: ', Nref[1], ', Nref NFE: ', Nref[2], ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 14)
p14
ggsave(file = paste0(dir_out, calc, '_gene_', file_out), plot = p14, height = 8, width = 16, units = 'in')



# t1e for new admixed AFR results-trying a Bar plot
# Doesn't work as well since you can't see the bar color when t1e is 0
p15 <- ggplot(results2 %>% filter(!(MACs == "Adjusted Ncc" | Method == "SKAT-O (External)" | Method == "SKAT-O (Internal + External)" |
                                      Method == "SKAT (External)" | Method == "SKAT (Internal + External)" |
                                      Method == "Burden (External)" | Method == "Burden (Internal + External)")),
              aes(x = Gene, y = Value, fill = Method, alpha = MACs)) +
  geom_bar(stat="identity", position="dodge", colour = "black") +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), width=.5, position=position_dodge(0.9)) +
  scale_fill_manual(values = colors3) +
  scale_alpha_manual(values=c(0.5, 1)) +
  facet_wrap(~Nref, ncol = 1, scales = 'fixed') +
  labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Internal Sample Size \nCases and Internal Controls: ', int_prune, '% pruned, 100% AFR \nCommon Controls: ', 
                                                ext_prune, '% pruned, 80% AFR, 20% NFE', 
                                                '\nNcc: ', Ncc, ', MAF: ', maf, ', Nref AFR: ', Nref[1], ', Nref NFE: ', Nref[2])) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 15)
p15
ggsave(file = paste0(dir_out, calc, '_BOXPLOT_gene_', file_out), plot = p15, height = 9, width = 18, units = 'in')

################################################################################
# t1e 100% NFE-proxECAT only
# t1e_int_v_ext = read.csv(paste0(dir, "T1e_NFE_99-80_maf", maf, ".csv"), header=T)
# t1e_ext_v_ext = read.csv(paste0(dir, "T1e_NFE_100v100-80v80_maf", maf, ".csv"), header=T)
# 
# t1e_int_v_ext = pivot_longer(t1e_int_v_ext, t1e_99:t1e_80, values_to="Value") %>%
#   mutate(Calculation = "Type I Error", MAF = maf, Pop = "100% NFE", Method = "ProxECAT",
#          Scenario = c("100% v 99%", "100% v 95%", "100% v 90%", "100% v 80%"),
#          Pruning = c("99%", "95%", "90%", "80%"))
# 
# t1e_ext_v_ext = pivot_longer(t1e_ext_v_ext, t1e_100v100:t1e_80v80, values_to="Value") %>%
#   mutate(Calculation = "Type I Error", MAF = maf, Pop = "100% NFE", Method = "ProxECAT",
#          Scenario = c("100% v 100%", "99% v 99%", "95% v 95%", "90% v 90%", "80% v 80%"),
#          Pruning = c("100%", "99%", "95%", "90%", "80%"))
# 
# results = rbind(t1e_int_v_ext, t1e_ext_v_ext)
# results$Scenario = factor(results$Scenario, levels=c("100% v 99%", "100% v 95%", "100% v 90%", 
#                                                      "100% v 80%", "100% v 100%",  "99% v 99%",
#                                                      "95% v 95%", "90% v 90%", "80% v 80%"))
# results$Pruning = factor(results$Pruning, levels=c("100%", "99%", "95%", "90%", "80%"))
# # get CI's
# results$Lower = '.'
# results$Upper = '.'
# # default level is 95% confidence
# nsim = 100
# for(i in 1:nrow(results)){
#   results$Lower[i] = binom.confint(nsim*results$Value[i], nsim, method=c("wilson"), type="central")$lower
#   results$Upper[i] = binom.confint(nsim*results$Value[i], nsim, method=c("wilson"), type="central")$upper
# }
# 
# results$Lower = as.numeric(results$Lower)
# results$Upper = as.numeric(results$Upper)
#############################################


# Proportion estimate plots
# t1e Scenario 1
# par(mfrow=c(2,3))
# plot(t1e.cc.prop.ests$maf_afr, t1e.case.prop.ests.s1$maf_afr, 
#      xlab = "AFR AF External Controls", ylab = "AFR AF Cases")
# abline(a=0, b=1)
# 
# plot(t1e.int.prop.ests$maf_afr, t1e.case.prop.ests.s1$maf_afr, 
#      xlab = "AFR AF Internal Controls", ylab = "AFR AF Cases")
# abline(a=0, b=1)
# 
# plot(t1e.cc.prop.ests$maf_afr, t1e.int.prop.ests.s1$maf_afr, 
#      xlab = "AFR AF External Controls", ylab = "AFR AF Internal Controls")
# abline(a=0, b=1)
# mtext("Scenario 1 Dataset Proportion Estimates-Type I error", side = 3, line = -3, outer = TRUE)


## Check to see if AFs differ greatly between AFR and NFE
# read in reference files
# hap.ref.afr.s1 = fread('chr19.block37.AFR.sim1.s1.ref.haps.gz')
# hap.ref.nfe.s1 = fread('chr19.block37.NFE.sim1.s1.ref.haps.gz')
# # Calculate AFs
# hap.ref.afr.s2 = fread('chr19.block37.AFR.sim1.s2.ref.haps.gz')
# hap.ref.nfe.s2 = fread('chr19.block37.NFE.sim1.s2.ref.haps.gz')
# count.ref.afr.s1 = data.frame(mac_afr=rowSums(hap.ref.afr.s1)) %>% mutate(maf_afr=mac_afr/(ncol(hap.ref.afr.s1)))
# count.ref.nfe.s1 = data.frame(mac_nfe=rowSums(hap.ref.nfe.s1)) %>% mutate(maf_nfe=mac_nfe/(ncol(hap.ref.nfe.s1)))
# count.ref.afr.s2 = data.frame(mac_afr=rowSums(hap.ref.afr.s2)) %>% mutate(maf_afr=mac_afr/(ncol(hap.ref.afr.s2)))
# count.ref.nfe.s2 = data.frame(mac_nfe=rowSums(hap.ref.nfe.s2)) %>% mutate(maf_nfe=mac_nfe/(ncol(hap.ref.nfe.s2)))
# # Plot resutls
# plot(count.ref.nfe.s1$maf_nfe, count.ref.afr.s1$maf_afr, xlim = c(0, 0.1), ylim = c(0, 0.1))
# plot(count.ref.nfe.s2$maf_nfe, count.ref.afr.s2$maf_afr, xlim = c(0, 0.1), ylim = c(0, 0.1))
