### Used to plot power results ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

# Pop = 'AFR'
# admx = '80-20'
Pop1 = 'AFR'
Pop2 = 'NFE'
admx_pop1 = 80
admx_pop2 = 20
Nsim = '42k' 
scen = 's2'
folder = '160v100v80'
pcase = 160
int_prune = 100
ext_prune = 80
Ncase = Nic = 5000
Ncc = 10000  
Nref_pop1 = 704
Nref_pop2 = 642
Neff = 10650
maf = 0.001 
sim_params = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_', Pop1, 'ref', Nref_pop1, '_', Pop2, 'ref', Nref_pop2)
sim_params1 = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_AFRref704_NFEref642')
sim_params2 = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_AFRref2000_NFEref2000')
sim_params3 = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_AFRref10000_NFEref10000')



pruning_plot = 'Separately and Sequentially-RAREsim v2.1.1' #Separately-RAREsim v2.1.1, Separately-R, Separately and Sequentially-RAREsim v2.1.1
# pruning = "pruneSepRaresim" #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
# data = 'by_gene'

# dir = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/Results/cc10k/'
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', Nsim, '_', Pop, '/', data, '/', folder, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', Pop1, '_', Pop2, '_pops/Sim_', Nsim, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', Pop1, '_', Pop2, '_pops/Sim_', Nsim, '/prox_gene_adj_', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params, '/prox_gene_adj_', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
dir1 = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params1, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
dir2 = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params2, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
dir3 = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params3, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')

# dir_out_t1e = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/typeI_error_plots/admixed/STATGEN/')
dir_out_power = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/power_plots/admixed/STATGEN/')

file_in = paste0(admx_pop1, Pop1, '_', admx_pop2, Pop2, '_', int_prune, 'v', ext_prune, '_maf', maf, '_', scen, '.csv')
file_out = paste0(admx_pop1, Pop1, "_", admx_pop2, Pop2, "_Nsim", Nsim, "_", folder, "_", int_prune, "v", ext_prune, "_Ncase", Ncase, "_Ncc", Ncc, "_maf", maf, "_", scen, '.jpg')


t1e_gene = read.csv(paste0(dir, "Power_all_gene_", file_in), header=T)

t1e_gene1 = read.csv(paste0(dir1, "Power_all_gene_", file_in), header=T)
t1e_gene2 = read.csv(paste0(dir2, "Power_all_gene_", file_in), header=T)
t1e_gene3 = read.csv(paste0(dir3, "Power_all_gene_", file_in), header=T)



t1e_gene1 = pivot_longer(t1e_gene1, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene2 = pivot_longer(t1e_gene2, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene3 = pivot_longer(t1e_gene3, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
# t1e_gene = pivot_longer(t1e_gene, prox_int:burden_all, names_to="Method", values_to="Value") %>% 
#   mutate(MAF = maf)

t1e_gene = rbind(t1e_gene1, t1e_gene2, t1e_gene3)

results2 = t1e_gene %>% mutate(MACs = rep(c("Unadjusted", "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Ncc ", Neff),
                                            "Unadjusted", "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Ncc ", Neff),
                                            "Unadjusted", "Unadjusted", "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Ncc ", Neff), paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Ncc ", Neff),
                                            "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Ncc ", Neff),
                                            "Unadjusted", "Unadjusted", "Unadjusted",
                                            "Unadjusted", "Unadjusted", "Unadjusted",
                                            "Unadjusted", "Unadjusted", "Unadjusted"), 12))

results2 = t1e_gene %>% mutate(MACs = rep(c("Unadjusted", "Unadjusted", "Adjusted", "Variant Adjusted Neff",
                                            "Unadjusted", "Unadjusted", "Adjusted", "Variant Adjusted Neff",
                                            "Unadjusted", "Unadjusted", "Unadjusted", "Adjusted", "Variant Adjusted Neff", "Adjusted", "Variant Adjusted Neff",
                                            "Unadjusted", "Adjusted", "Variant Adjusted Neff",
                                            "Unadjusted", "Unadjusted", "Unadjusted",
                                            "Unadjusted", "Unadjusted", "Unadjusted",
                                            "Unadjusted", "Unadjusted", "Unadjusted"), 36))

# results2 = t1e_gene %>% mutate(MACs = rep(c("Unadjusted", "Unadjusted", "Variant Adjusted Ncc", "Variant Adjusted Neff", 
#                                             "Gene Adjusted Ncc", "Gene Adjusted Neff"), 72))

results2$Method = factor(results2$Method, levels=c("prox_int", "prox_ext", "prox_ext_adj_Ncc", "prox_ext_adj_Neff", 
                                                   "proxW_int", "proxW_ext", "proxW_ext_adj_Ncc", "proxW_ext_adj_Neff",
                                                   "prox2_int", "prox2_ext", "prox2_all", "prox2_ext_adj_Ncc", "prox2_ext_adj_Neff", "prox2_all_adj_Ncc", "prox2_all_adj_Neff",
                                                   "iecat_all", "iecat_all_adj_Ncc", "iecat_all_adj_Neff",
                                                   "skato_int", "skato_ext", "skato_all",
                                                   "skat_int", "skat_ext", "skat_all",
                                                   "burden_int", "burden_ext", "burden_all"),
                         labels=rep(c("ProxECAT", "ProxECAT-weighted", "LogProx", "iECAT-O", "SKAT-O", "SKAT", "Burden"), times=c(4, 4, 7, 3, 3, 3, 3)))

# results2$Method = factor(results2$Method, levels=c("prox_int", "prox_ext", "prox_ext_var_adj_Ncc", "prox_ext_var_adj_Neff", "prox_ext_gene_adj_Ncc", "prox_ext_gene_adj_Neff", 
#                                                    "proxW_int", "proxW_ext", "proxW_ext_var_adj_Ncc", "proxW_ext_var_adj_Neff", "proxW_ext_gene_adj_Ncc", "proxW_ext_gene_adj_Neff"),
#                          labels=rep(c("ProxECAT", "ProxECAT-weighted"), times=c(6, 6)))

results2 = results2 %>% mutate(Data = rep(c("Internal", "External", "External", "External",
                                            "Internal", "External", "External", "External",
                                            "Internal", "External", "Internal + External", "External", "External", "Internal + External", "Internal + External",
                                            "Internal + External", "Internal + External", "Internal + External",
                                            "Internal", "External", "Internal + External",
                                            "Internal", "External", "Internal + External",
                                            "Internal", "External", "Internal + External"), 36))

# results2 = results2 %>% mutate(Data = rep(c("Internal", "External", "External", "External", "External", "External"), 72))
results2 = results2 %>% mutate(Nref = rep(c("AFR: 704, NFE: 642", "AFR: 2000, NFE: 2000", "AFR: 10000, NFE: 10000"), each=324))

# results2$MAF = factor(results2$MAF, levels = c(0.01, 0.001))
# results2$Calculation = factor(results2$Calculation, levels = c("Type I Error", "Power"))
# results2$Scenario = factor(results2$Scenario)
results2$MAF = factor(results2$MAF)
results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Ncc ", Neff)))
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Adjusted", "Variant Adjusted Neff"))
results2$Gene = factor(results2$Gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
                                               "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))
results2$Nref = factor(results2$Nref, levels=c("AFR: 704, NFE: 642", "AFR: 2000, NFE: 2000", "AFR: 10000, NFE: 10000"))

# results2$MAF = factor(results2$MAF)
# results2$Data = factor(results2$Data, levels=c("Internal", "External"))
# results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Variant Adjusted Ncc", "Variant Adjusted Neff", "Gene Adjusted Ncc", "Gene Adjusted Neff"))
# results2$Gene = factor(results2$Gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
#                                                "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))
# results2$Nref = factor(results2$Nref, levels=c(500, 2000, 10000))

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
colors3 = c("#009E73", "#0072B2", "#D55E00")
colors_adj = c("#009E73", "#0072B2", "#56B4E9", "#D55E00", "#E69F00")
colors_nfe = c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999",
               "lightgreen", "#56B4E9", "#E69F00", "pink")
colors_all_nfe = "#009E73"
colors_meth = c("#56B4E9", "#0072B2", "#009E73", "#CC79A7", "#BC9F4C", "#E69F00", "#D55E00")
colors_meth_adj = c("#56B4E9", "#0072B2", "#009E73", "#CC79A7")
colors_gene = c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999",
                "lightgreen", "#56B4E9", "#E69F00", "#BC9F4C", "pink", "red", "purple")
colors_prox = c("#009E73", "#0072B2", "#56B4E9")

# Controls everything in the graph
# base_size = 25


# t1e for admixed by gene results from pcasev100vpconf pipeline-a for all methods
p10a <- ggplot(results2, aes(x=Gene, y=Value, color=Method, shape=MACs)) +
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
                                                '\nCases and Internal Controls: 100% AFR, External Controls: 80% AFR 20% NFE \nNsim: ', Nsim, ', Ncase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', Nref ', Pop1, ': ', Nref_pop1, ', Nref ', Pop2, ': ', Nref_pop2, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 15)
p10a
ggsave(file = paste0(dir_out_t1e, file_out), plot = p10a, height = 8, width = 16, units = 'in')

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
  scale_color_manual(values=c("#009E73", "#CC79A7")) +
  scale_shape_manual(values = c(1, 16)) +
  # scale_linetype_manual(values = c("solid", "dotted", "dashed")) +
  # scale_size_manual(values = c(0.5, 1, 3)) +
  facet_wrap(~Nref, ncol = 1, scales = 'fixed') +
  # facet_wrap(~Data, ncol = 1) +
  labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene and Reference Sample Size: Cases vs Internal + Common Controls \nCases and Internal Controls: ', int_prune, '% pruned, 100% AFR \nCommon Controls: ', 
                                                ext_prune, '% pruned, 80% AFR, 20% NFE', 
                                                '\nNsim: ', Nsim, ', Ncase: ', Ncase, ', Nic: ', Nic, ', Ncc: ', Ncc, ', MAF: ', maf)) +
  # '\nPop: Admixed ', admx, " ", Pop1, '-', Pop2, ', Nsim: ', Nsim, ', Ncase: ', Ncase, ', Ncc: ', Ncc, ', Nref: ', Nref, ', MAF: 0.001')) +
  # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
  theme_bw(base_size = 15)
p13
ggsave(file = paste0(dir_out_t1e, 't1e_gene_Nref_int_and_ext_only_', file_out), plot = p13, height = 8, width = 16, units = 'in')
