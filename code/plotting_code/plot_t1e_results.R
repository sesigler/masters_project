### Used to plot type I error results ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

calc = 'T1e'
Pops = c('AFR', 'NFE')
admx_props = c(80, 20)
scen = 's2'
sub_scen = 'default'
int_prune = 100
ext_prune = 80
Ncase = Nic = 2000
Ncc= 10000
Nref = 2000
maf = 0.001 
Ncc = 10000
Neff = 10046 #median value


dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/')
dir_out = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/typeI_error_plots/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/')
# dir_out_power = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/power_plots/', data, '/')

file_in = paste0(scen, "_", sub_scen, "_maf", maf, ".csv")
file_out = paste0(scen, "_", sub_scen, "_maf", maf, '.jpg')

# read in the results
# Proportion Estimates 
# t1e_case_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_case_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_cc_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_cc_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_int_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_int_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

results = read.csv(paste0(dir, calc, "_all_gene_", file_in), header=T)
# t1e_gene = read.csv(paste0(dir, "T1e_all_prox_gene_adj_", file_in), header=T)
t1e_gene1 = read.csv(paste0(dir1, "T1e_all_gene_", file_in), header=T)
t1e_gene2 = read.csv(paste0(dir2, "T1e_all_gene_", file_in), header=T)
t1e_gene3 = read.csv(paste0(dir3, "T1e_all_gene_", file_in), header=T)


# puts in a format for ggplot
t1e_gene1 = pivot_longer(t1e_gene1, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene2 = pivot_longer(t1e_gene2, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene3 = pivot_longer(t1e_gene3, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)
t1e_gene = rbind(t1e_gene1, t1e_gene2, t1e_gene3)

results = pivot_longer(results, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(MAF = maf)



results2 = results %>% mutate(MACs = rep(c("Unadjusted", "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Neff ", Neff),
                                            "Unadjusted", "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Neff ", Neff),
                                            "Unadjusted", "Unadjusted", "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Neff ", Neff), paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Neff ", Neff),
                                            "Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Neff ", Neff),
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
                                            "Internal", "External", "Internal + External"), 12))

# results2 = results2 %>% mutate(Data = rep(c("Internal", "External", "External", "External", "External", "External"), 72))
results2 = results2 %>% mutate(Nref = rep(c("AFR: 704, NFE: 642", "AFR: 2000, NFE: 2000", "AFR: 10000, NFE: 10000"), each=324))

# results2$MAF = factor(results2$MAF, levels = c(0.01, 0.001))
# results2$Calculation = factor(results2$Calculation, levels = c("Type I Error", "Power"))
# results2$Scenario = factor(results2$Scenario)
results2$MAF = factor(results2$MAF)
results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), paste0("Variant Adjusted Neff ", Neff)), 
                       labels = c("Unadjusted", paste0("Variant Adjusted Ncc ", Ncc), "Adjusted"))
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

# Type I error graph
p1 <- ggplot(results2, aes(x=Method, y=Value, color=Pop)) +
        geom_point(aes(shape=Pop), size=5, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
        geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors_all_nfe) +
        facet_grid(~Data, scales="free", space="free") +
        labs(y='Type I Error', x='Method', title=paste0('Scenario 1: Type I Error 100% vs ', ext_prune, '% (10k cc) \nMAF=0.001')) +
        theme_bw(base_size = 20)
p1
ggsave(file = paste0(dir_out, 't1e_', Pop2, '_', int_prune, '_v_', ext_prune, '.jpg'),
       plot = p1, height = 8, width = 15, units = 'in')

#Pruning NFE plot
p2 <- ggplot(results, aes(x=Pruning, y=Value, color=Scenario)) +
        geom_point(shape=rep(c(16, 17), times=c(4, 5)), size=5, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
        geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors_nfe) +
        # facet_grid(~Data, scales="free", space="free") +
        labs(y='Type I Error', x='Pruning', title='Type I Error for Various Pruning Scenarios of Internal Cases \nand External Controls \nPop=100% NFE, MAF=0.001') +
        guides(color = guide_legend(override.aes=list(shape = rep(c(16, 17), times=c(4, 5))))) +
        # guides(fill=guide_legend(title="Pruning Scenario for Internal Cases vs External Controls")) +
        theme_bw(base_size = 20)
p2
ggsave(file = paste0(dir_out, 't1e_', Pop2, '_all_pruning_scenarios.jpg'),
       plot = p2, height = 8, width = 12, units = 'in')

# Fun or syn variants used all methods
p3 <- ggplot(results2, aes(x=Method, y=Value, color=Variants_Used)) +
        geom_point(aes(shape=Variants_Used), size=5, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
        geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors3) +
        facet_grid(~Data, scales="free", space="free") +
        labs(y='Type I Error', x='Method', title=paste0('Type I Error: 100% vs ', ext_prune, '% (10k cc) \nPop=100% NFE, MAF=0.001')) +
        theme_bw(base_size = 20)
p3
ggsave(file = paste0(dir_out, 't1e_', Pop2, '_', int_prune, '_v_', ext_prune, '_variants_used.jpg'),
       plot = p3, height = 8, width = 15, units = 'in')

# t1e for OG hap pruning straight down to 100% fun and 100% syn
p4 <- ggplot(results2, aes(x=Method, y=Value, color=Variants_Used)) +
        geom_point(aes(shape=Variants_Used), size=5, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
        geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors3) +
        facet_grid(~Data, scales="free", space="free") +
        labs(y='Type I Error', x='Method', title=paste0('Type I Error: 100% vs ', ext_prune, '% (10k cc) \nOriginal Hapgen Haplotype Pruned Down\nPop=100% NFE, MAF=0.001')) +
        theme_bw(base_size = 20)
p4
ggsave(file = paste0(dir_out, 't1e_og_hap_', Pop2, '_', int_prune, '_v_', ext_prune, '_variants_used.jpg'),
       plot = p4, height = 8, width = 15, units = 'in')

# t1e for different pruning scenarios for 
p5 <- ggplot(results2, aes(x=Method, y=Value, color=Variants_Used)) +
        geom_point(aes(shape=Variants_Used), size=5, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
        geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors3) +
        facet_grid(~Data, scales="free", space="free") +
        labs(y='Type I Error', x='Method', title=paste0('Type I Error: ', int_prune, '% vs ', ext_prune, '% (10k cc) \nPruning: ', pruning, ' (RAREsim v2.1.1-100% fun, R-100% syn)\nPop=100% NFE, MAF=0.001')) +
        theme_bw(base_size = 20)
p5
ggsave(file = paste0(dir_out, 't1e_', Pop2, '_', pruning, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
       plot = p5, height = 8, width = 15, units = 'in')

# t1e for different pruning scenarios
p6 <- ggplot(results2, aes(x=Method, y=Value, color=Variants_Used)) +
        geom_point(aes(shape=Variants_Used), size=5, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
        geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors3) +
        facet_grid(~Data, scales="free", space="free") +
        labs(y='Type I Error', x='Method', title=paste0('Type I Error: ', int_prune, '% vs ', ext_prune, '% from ', folder, ' Data (10k cc) \nPruning: ', pruning_plot, '\nPop=100% NFE, MAF=0.001')) +
        # labs(y='Type I Error', x='Method', title=paste0('Type I Error: ', int_prune, '% vs ', ext_prune, '% from ', data, ' Data (10k cc) \nPruning: ', pruning_plot, '\nPop=100% NFE, MAF=0.001')) +
        theme_bw(base_size = 15)
p6
ggsave(file = paste0(dir_out, 't1e_all_', Pop2, '_', pruning, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
       plot = p6, height = 8, width = 16, units = 'in')
# ggsave(file = paste0(dir_out, 't1e_all_', Pop2, '_', pruning, '_', data, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
#        plot = p6, height = 8, width = 16, units = 'in') #int v ext

# t1e for by gene results
p7 <- ggplot(results2, aes(x=Gene, y=Value, color=Method)) +
        geom_point(aes(shape=Method), size=3, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
        # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        # scale_y_continuous(breaks=c(0, 0.05, 0.20, 0.40, 0.60)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors_meth) +
        scale_shape_manual(values = c(16, 18, 17, 15, 9, 10, 12)) +
        facet_wrap(~Data, ncol = 1, scales = 'free_y') +
        labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene: ', int_prune, '% vs ', ext_prune, '% from ', folder, ' Data (10k cc) \nPruning: ', pruning_plot, '\nPop=100% NFE, MAF=0.001')) +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
        theme_bw(base_size = 15)
p7
ggsave(file = paste0(dir_out, 't1e_gene_', Pop2, '_', pruning, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
       plot = p7, height = 8, width = 16, units = 'in')

# t1e for by gene results from 120v100v80 pipeline
p8 <- ggplot(results2 %>% filter(Calculation == "Type I Error"), aes(x=Gene, y=Value, color=Method)) +
        geom_point(aes(shape=Method), size=3, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.05, linetype=2, linewidth=1) +
        # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        # scale_y_continuous(breaks=c(0, 0.05, 0.25, 0.5, 0.75, 1)) +
        # scale_y_continuous(breaks=c(0, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors_meth) +
        scale_shape_manual(values = c(16, 18, 17, 15, 9, 10, 12)) +
        facet_wrap(~Data, ncol = 1, scales = 'free_y') +
        # facet_wrap(~Data, ncol = 1) +
        labs(y='Type I Error', x='Gene', title=paste0('Type I Error by Gene: ', int_prune, '% vs ', ext_prune, '% from ', folder, ' Data (10k cc) \nPruning: ', pruning_plot, '\nPop=100% ', Pop, ', MAF=0.001')) +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
        theme_bw(base_size = 15)
p8
ggsave(file = paste0(dir_out_t1e, 't1e_gene_', Nsim, '_', Pop, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
       plot = p8, height = 8, width = 16, units = 'in')
# ggsave(file = paste0(dir_out_t1e, 't1e_gene_same_y_', Pop2, '_', pruning, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
#        plot = p8, height = 8, width = 16, units = 'in')

# power for by gene results from 120v100v80 pipeline
p9 <- ggplot(results2 %>% filter(Calculation == "Power") %>%
               filter(!(Method == "SKAT-O" & (Data == "External" | Data == "Internal + External"))) %>%
               filter(!(Method == "SKAT" & (Data == "External" | Data == "Internal + External"))) %>%
               filter(!(Method == "Burden" & (Data == "External" | Data == "Internal + External"))), 
             aes(x=Gene, y=Value, color=Method)) +
        geom_point(aes(shape=Method), size=3, position=position_dodge(width=0.5)) +
        geom_hline(yintercept=0.80, linetype=2, linewidth=1) +
        # geom_hline(yintercept=1, linetype="blank", linewidth=1.5) +
        # scale_y_continuous(limits=c(0, 1)) +
        scale_y_continuous(breaks=c(0, 0.20, 0.40, 0.60, 0.80, 1)) +
        # scale_y_continuous(breaks=c(0, 0.05, 0.20, 0.40, 0.60)) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1, width=.3, position=position_dodge(width=0.5)) +
        scale_color_manual(values=colors_meth) +
        scale_shape_manual(values = c(16, 18, 17, 15, 9, 10, 12)) +
        # facet_wrap(~Data, ncol = 1, scales = 'free_y') +
        facet_wrap(~Data, ncol = 1) +
        labs(y='Power', x='Gene', title=paste0('Power by Gene: ', pcase, '% Cases vs 100% Internal Controls vs ', ext_prune, '% External Controls from ', folder, ' Data (10k cc) \nPruning: ', pruning_plot, '\nPop=100% ', Pop, ', MAF=0.001')) +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
        theme_bw(base_size = 15)
p9
ggsave(file = paste0(dir_out_power, 'power_gene_', Nsim, '_', Pop, '_', folder, '_maf', maf, '.jpg'),
       plot = p9, height = 8, width = 15, units = 'in')
 

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
