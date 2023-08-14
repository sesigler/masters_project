### Used to plot histograms of p-values ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread
library(patchwork)

Pop1 = "AFR"
Pop2 = "NFE"
scen = "s1"
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc10k'  #Number of common controls: 'cc5k' or 'cc10k'
int_prune = 100
ext_prune = 100

dir = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/Results/cc10k/'

### Type 1 error
t1e_homo = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e_adj = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

# puts in a format for ggplot
t1e_homo = pivot_longer(t1e_homo, prox_p:skat_all_p, names_to="Method", values_to="P-Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, 
         Pop = "Homogeneous", Configuration = "Homogeneous-Unadjusted",
         Data  = rep(c("External", "Internal", "External", "Internal + External", "Internal",
                       "Internal + External", "Internal", "External", "Internal + External"), 100))

t1e = pivot_longer(t1e, prox_p:skat_all_p, names_to="Method", values_to="P-Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, 
         Pop = "Admixed", Configuration = "Admixed-Unadjusted",
         Data  = rep(c("External", "Internal", "External", "Internal + External", "Internal",
                       "Internal + External", "Internal", "External", "Internal + External"), 100))

t1e_adj = pivot_longer(t1e_adj, prox_p_adj:iecat_p_adj, names_to="Method", values_to="P-Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, 
         Pop = "Admixed", Configuration = "Admixed-Adjusted",
         Data  = rep(c("External", "External", "Internal + External", "Internal + External"), 100))

results = rbind(t1e_homo, t1e, t1e_adj)

results$Method = factor(results$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
                                                   "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p"),
                         labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
                                  "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O"))

results$Calculation = factor(results$Calculation)
results$Scenario = factor(results$Scenario)
results$MAF = factor(results$MAF)
results$Pop = factor(results$Pop, levels=c("Admixed", "Homogeneous"))
results$Data = factor(results$Data, levels=c("Internal", "External", "Internal + External"))
results$Configuration = factor(results$Configuration, levels=c("Homogeneous-Unadjusted", "Admixed-Unadjusted", "Admixed-Adjusted"))

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors = c("#009E73", "#0072B2", "#D55E00")

# P-value Distribution
ggplot(results,
       aes(x=`P-Value`, color=Configuration, fill=Configuration)) +
  geom_histogram(alpha=0.5, binwidth = 0.05) +
  #geom_histogram(aes(y=after_stat(density)), alpha=0.5, binwidth = 0.05) +
  #geom_density(alpha=0.2) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  facet_grid(Data~Method, space = "free_x", scales = "free_x")

hist1 <- results %>%
  filter(Method != "iECAT-O") %>%
  ggplot(results, aes(x=`P-Value`, color=Configuration, fill=Configuration)) +
    geom_histogram(alpha=0.5, binwidth = 0.05) +
    #geom_histogram(aes(y=after_stat(density)), alpha=0.5, binwidth = 0.05) +
    #geom_density(alpha=0.2) +
    scale_color_manual(values=colors)+
    scale_fill_manual(values=colors) +
    facet_grid(rows = vars(Data), cols = vars(Method))

ggplot(results %>% filter(Method == "iECAT-O"),
       aes(x=`P-Value`, color=Configuration, fill=Configuration)) +
  geom_histogram(alpha=0.5, binwidth = 0.05) +
  #geom_histogram(aes(y=after_stat(density)), alpha=0.5, binwidth = 0.05) +
  #geom_density(alpha=0.2) +
  scale_color_manual(values=colors)+
  scale_fill_manual(values=colors) +
  facet_grid(Data~Method)

