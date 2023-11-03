### Used to plot type I error results ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

Pop1 = 'AFR'
Pop2 = 'NFE'
# scen = 's1'
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc10k'  #Number of common controls: 'cc5k' or 'cc10k'
int_prune = 100
ext_prune = 99
pruning = "pruneSepRaresim" #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '100v99'
pruning_plot = 'Separately-RAREsim v2.1.1' #Separately-RAREsim v2.1.1, Separately-R
data = 'by_gene'

# dir = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/Results/cc10k/'
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/', int_prune, 'v', ext_prune, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', folder, '/', int_prune, 'v', ext_prune, '/')
dir_out = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Results/typeI_error_plots/', pruning, '/', data, '/')

# read in the results
# Proportion Estimates 
# t1e_case_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_case_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_cc_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_cc_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_int_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_int_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

# t1e
# t1e_all = read.csv(paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# t1e_all_adj = read.csv(paste0(dir, "T1e_all_adj_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# t1e_all_homo = read.csv(paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# t1e_all_pruning = read.csv(paste0(dir, "T1e_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# t1e_all_pruning = read.csv(paste0(dir, "T1e_all_", pruning, "_", data, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# t1e_skat_pruning = read.csv(paste0(dir, "T1e_all_skat_", pruning, "_", data, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# t1e_skat_pruning = read.csv(paste0(dir, "T1e_all_skat_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
t1e_gene = read.csv(paste0(dir, "T1e_all_gene_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)

# puts in a format for ggplot
# t1e_all = pivot_longer(t1e_all, prox_p:skat_all_p, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "Admixed")
# t1e_all_adj = pivot_longer(t1e_all_adj, prox_p_adj:iecat_p_adj, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "Admixed")
# t1e_all_homo = pivot_longer(t1e_all_homo, prox_p:skat_all_p, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "100% NFE")

# t1e_all_skat = pivot_longer(t1e_all_skat, iecat_p_syn:skat_all_p_syn, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "100% NFE")
# t1e_all_og_hap = pivot_longer(t1e_all_og_hap, prox_p:skat_all_p_syn, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "100% NFE")

# t1e_all_pruning = pivot_longer(t1e_all_pruning, prox_p:skat_all_p_syn, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "100% NFE")
# t1e_skat_pruning = pivot_longer(t1e_skat_pruning, skat_int_p:burden_all_p_syn, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "100% NFE")
# t1e_all_pruning = pivot_longer(t1e_all_pruning, prox_p:burden_all_p_syn, names_to="Method", values_to="Value") %>%
#   mutate(Calculation = "Type I Error", MAF = maf, Pop = "100% NFE")
t1e_gene = pivot_longer(t1e_gene, prox_int:burden_all, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", MAF = maf, Pop = "100% NFE")


# results = rbind(t1e_all_homo, t1e_all, t1e_all_adj)
# results = t1e_all_homo
# results = rbind(t1e_all_homo, t1e_all_skat)
# results = t1e_all_pruning
# results = rbind(t1e_all_pruning, t1e_skat_pruning)
# results_all = t1e_all_pruning
# results_skat = t1e_skat_pruning
results = t1e_gene

# results2 = results %>% mutate(Data = c("External", "Internal", "External", "Internal + External", "Internal", 
#                                        "Internal + External", "Internal", "External", "Internal + External",
#                                        "External", "Internal", "External", "Internal + External", "Internal", 
#                                        "Internal + External", "Internal", "External", "Internal + External",
#                                        "External", "External", "Internal + External", "Internal + External")) %>%
#   mutate(MACs = rep(c("Unadjusted", "Adjusted"), c(18, 4))) %>%
#   mutate(Configuration = rep(c("Homogeneous-Unadjusted", "Admixed-Unadjusted", "Admixed-Adjusted"), c(9, 9, 4)))

# results2 = results %>% mutate(Data = c("External", "Internal", "External", "Internal + External", "Internal", 
#                                        "Internal + External", "Internal", "External", "Internal + External"))

# results2 = results %>% mutate(Data = c("External", "Internal", "External", "Internal + External", "Internal", 
#                                        "Internal + External", "Internal", "External", "Internal + External",
#                                        "Internal + External", "Internal", "External", "Internal + External")) %>% 
#   mutate(Variants_Used = rep(c("Functional and Synonymous", "Functional Only", "Synonymous Only"), times=c(5, 4, 4)))

# results_all = results_all %>% mutate(Data = c("External", "Internal", "External", "Internal + External", "Internal", 
#                                               "Internal + External", "Internal", "External", "Internal + External",
#                                               "Internal + External", "Internal", "External", "Internal + External")) %>% 
#   mutate(Variants_Used = rep(c("Functional and Synonymous", "Functional Only", "Synonymous Only"), times=c(5, 4, 4)))
# 
# results_skat = results_skat%>% mutate(Data = rep(c("Internal", "External", "Internal + External"), 4)) %>% 
#   mutate(Variants_Used = rep(c("Functional Only", "Synonymous Only"), times=c(6, 6)))

# results2 = results %>% mutate(Data = c("External", "Internal", "External", "Internal + External", "Internal", 
#                                        "Internal + External", "Internal", "External", "Internal + External",
#                                        "Internal + External", "Internal", "External", "Internal + External",
#                                        "Internal", "External", "Internal + External",
#                                        "Internal", "External", "Internal + External",
#                                        "Internal", "External", "Internal + External",
#                                        "Internal", "External", "Internal + External")) %>%
#   mutate(Variants_Used = rep(c("Functional and Synonymous", "Functional Only", "Synonymous Only", 
#                                "Functional Only", "Synonymous Only"), times=c(5, 4, 4, 6, 6)))

# results2$Method = factor(results2$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
#                                                    "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p",
#                                                    "prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
#                                                    "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p",
#                                                    "prox_p_adj", "prox2_p_adj", "prox2_all_p_adj", "iecat_p_adj"),
#                          labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
#                                   "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O",
#                                   "ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
#                                   "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O",
#                                   "ProxECAT", "LogProx", "LogProx", "iECAT-O"))

# results2$Method = factor(results2$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
#                                                    "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p"),
#                          labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
#                                   "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O"))

# results2$Method = factor(results2$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
#                                                    "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p",
#                                                    "iecat_p_syn", "skat_int_p_syn", "skat_ext_p_syn", "skat_all_p_syn"),
#                          labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
#                                   "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O",
#                                   "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O"))

# results_all$Method = factor(results_all$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
#                                                          "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p",
#                                                          "iecat_p_syn", "skat_int_p_syn", "skat_ext_p_syn", "skat_all_p_syn"),
#                             labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
#                                      "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O",
#                                      "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O"))
# 
# results_skat$Method = factor(results_skat$Method, levels = c("skat_int_p", "skat_ext_p", "skat_all_p",
#                                                              "burden_int_p", "burden_ext_p", "burden_all_p",
#                                                              "skat_int_p_syn", "skat_ext_p_syn", "skat_all_p_syn",
#                                                              "burden_int_p_syn", "burden_ext_p_syn", "burden_all_p_syn"),
#                              labels=rep(c("SKAT", "SKAT", "SKAT", "Burden", "Burden", "Burden"), 2))

# results2$Method = factor(results2$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
#                                                    "iecat_p", "skato_int_p", "skato_ext_p", "skato_all_p",
#                                                    "iecat_p_syn", "skato_int_p_syn", "skato_ext_p_syn", "skato_all_p_syn",
#                                                    "skat_int_p", "skat_ext_p", "skat_all_p",
#                                                    "burden_int_p", "burden_ext_p", "burden_all_p",
#                                                    "skat_int_p_syn", "skat_ext_p_syn", "skat_all_p_syn",
#                                                    "burden_int_p_syn", "burden_ext_p_syn", "burden_all_p_syn"),
#                          labels=rep(c("proxECAT", "LogProx", "iECAT-O", "SKAT-O", "iECAT-O", "SKAT-O",
#                                       "SKAT", "Burden", "SKAT", "Burden"), times=c(2, 3, 1, 3, 1, 3, 3, 3, 3, 3)))

results2 = results %>% mutate(Data = rep(c("Internal", "External", "Internal", "External",
                                           "Internal", "External", "Internal + External",
                                           "Internal + External", "Internal", "External", "Internal + External",
                                           "Internal", "External", "Internal + External",
                                           "Internal", "External", "Internal + External"), 12))

results2$Method = factor(results2$Method, levels=c("prox_int", "prox_ext", "proxW_int", "proxW_ext",
                                                   "prox2_int", "prox2_ext", "prox2_all",
                                                   "iecat_all", "skato_int", "skato_ext", "skato_all",
                                                   "skat_int", "skat_ext", "skat_all",
                                                   "burden_int", "burden_ext", "burden_all"),
                         labels=rep(c("ProxECAT", "ProxECAT-weighted", "LogProx", "iECAT-O", "SKAT-O", 
                                      "SKAT", "Burden"), times=c(2, 2, 3, 1, 3, 3, 3)))



# results2 = rbind(results_all, results_skat)

#results2$MAF = factor(results2$MAF, levels = c(0.01, 0.001))
results2$Calculation = factor(results2$Calculation)
# results2$Scenario = factor(results2$Scenario)
results2$MAF = factor(results2$MAF)
# results2$Pop = factor(results2$Pop, levels=c("Admixed", "Homogeneous"))
results2$Pop = factor(results2$Pop)
results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))
# results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Adjusted"))
# results2$Configuration = factor(results2$Configuration, levels=c("Homogeneous-Unadjusted", "Admixed-Unadjusted", "Admixed-Adjusted"))
# results2$Variants_Used = factor(results2$Variants_Used, levels=c("Functional and Synonymous", "Functional Only", "Synonymous Only"))
results2$Gene = factor(results2$Gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
                                               "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))
# results2$Method = factor(results2$Method, levels = c("ProxECAT", "ProxECAT-weighted", "LogProx", "iECAT-O",
#                                                      "SKAT-O", "SKAT", "Burden"))

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
colors_nfe = c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999",
               "lightgreen", "#56B4E9", "#E69F00", "pink")
colors_all_nfe = "#009E73"
colors_meth = c("#56B4E9", "#0072B2", "#009E73", "#CC79A7", "#BC9F4C", "#E69F00", "#D55E00")
colors_gene = c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#999999",
                "lightgreen", "#56B4E9", "#E69F00", "#BC9F4C", "pink", "red", "purple")

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

# t1e for by gene resutls
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
        labs(y='Type I Error', x='Method', title=paste0('Type I Error by Gene: ', int_prune, '% vs ', ext_prune, '% from ', folder, ' Data (10k cc) \nPruning: ', pruning_plot, '\nPop=100% NFE, MAF=0.001')) +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65))
        theme_bw(base_size = 15)
p7
ggsave(file = paste0(dir_out, 't1e_gene_', Pop2, '_', pruning, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'),
       plot = p7, height = 8, width = 16, units = 'in')



################################################################################
# t1e 100% NFE-proxECAT only
t1e_int_v_ext = read.csv(paste0(dir, "T1e_NFE_99-80_maf", maf, ".csv"), header=T)
t1e_ext_v_ext = read.csv(paste0(dir, "T1e_NFE_100v100-80v80_maf", maf, ".csv"), header=T)

t1e_int_v_ext = pivot_longer(t1e_int_v_ext, t1e_99:t1e_80, values_to="Value") %>%
  mutate(Calculation = "Type I Error", MAF = maf, Pop = "100% NFE", Method = "ProxECAT",
         Scenario = c("100% v 99%", "100% v 95%", "100% v 90%", "100% v 80%"),
         Pruning = c("99%", "95%", "90%", "80%"))

t1e_ext_v_ext = pivot_longer(t1e_ext_v_ext, t1e_100v100:t1e_80v80, values_to="Value") %>%
  mutate(Calculation = "Type I Error", MAF = maf, Pop = "100% NFE", Method = "ProxECAT",
         Scenario = c("100% v 100%", "99% v 99%", "95% v 95%", "90% v 90%", "80% v 80%"),
         Pruning = c("100%", "99%", "95%", "90%", "80%"))

results = rbind(t1e_int_v_ext, t1e_ext_v_ext)
results$Scenario = factor(results$Scenario, levels=c("100% v 99%", "100% v 95%", "100% v 90%", 
                                                     "100% v 80%", "100% v 100%",  "99% v 99%",
                                                     "95% v 95%", "90% v 90%", "80% v 80%"))
results$Pruning = factor(results$Pruning, levels=c("100%", "99%", "95%", "90%", "80%"))
# get CI's
results$Lower = '.'
results$Upper = '.'
# default level is 95% confidence
nsim = 100
for(i in 1:nrow(results)){
  results$Lower[i] = binom.confint(nsim*results$Value[i], nsim, method=c("wilson"), type="central")$lower
  results$Upper[i] = binom.confint(nsim*results$Value[i], nsim, method=c("wilson"), type="central")$upper
}

results$Lower = as.numeric(results$Lower)
results$Upper = as.numeric(results$Upper)
#############################################


# Proportion estimate plots
# t1e Scenario 1
par(mfrow=c(2,3))
plot(t1e.cc.prop.ests$maf_afr, t1e.case.prop.ests.s1$maf_afr, 
     xlab = "AFR AF External Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(t1e.int.prop.ests$maf_afr, t1e.case.prop.ests.s1$maf_afr, 
     xlab = "AFR AF Internal Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(t1e.cc.prop.ests$maf_afr, t1e.int.prop.ests.s1$maf_afr, 
     xlab = "AFR AF External Controls", ylab = "AFR AF Internal Controls")
abline(a=0, b=1)
mtext("Scenario 1 Dataset Proportion Estimates-Type I error", side = 3, line = -3, outer = TRUE)


## Check to see if AFs differ greatly between AFR and NFE
# read in reference files
hap.ref.afr.s1 = fread('chr19.block37.AFR.sim1.s1.ref.haps.gz')
hap.ref.nfe.s1 = fread('chr19.block37.NFE.sim1.s1.ref.haps.gz')
# Calculate AFs
hap.ref.afr.s2 = fread('chr19.block37.AFR.sim1.s2.ref.haps.gz')
hap.ref.nfe.s2 = fread('chr19.block37.NFE.sim1.s2.ref.haps.gz')
count.ref.afr.s1 = data.frame(mac_afr=rowSums(hap.ref.afr.s1)) %>% mutate(maf_afr=mac_afr/(ncol(hap.ref.afr.s1)))
count.ref.nfe.s1 = data.frame(mac_nfe=rowSums(hap.ref.nfe.s1)) %>% mutate(maf_nfe=mac_nfe/(ncol(hap.ref.nfe.s1)))
count.ref.afr.s2 = data.frame(mac_afr=rowSums(hap.ref.afr.s2)) %>% mutate(maf_afr=mac_afr/(ncol(hap.ref.afr.s2)))
count.ref.nfe.s2 = data.frame(mac_nfe=rowSums(hap.ref.nfe.s2)) %>% mutate(maf_nfe=mac_nfe/(ncol(hap.ref.nfe.s2)))
# Plot resutls
plot(count.ref.nfe.s1$maf_nfe, count.ref.afr.s1$maf_afr, xlim = c(0, 0.1), ylim = c(0, 0.1))
plot(count.ref.nfe.s2$maf_nfe, count.ref.afr.s2$maf_afr, xlim = c(0, 0.1), ylim = c(0, 0.1))
