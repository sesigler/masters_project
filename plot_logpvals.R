# Plot the -log10pvalues for two different pruning scenarios

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

Pop1 = "AFR"
Pop2 = "NFE"
scen = "s1"
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
int_prune = 100
ext_prune = 99
ext_prune2 = 80

dir = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/Results/cc10k/'

# Read in p-values
# t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", int_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_conf = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e_conf = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune2, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)


# puts in a format for ggplot
t1e = pivot_longer(t1e, prox_p:skat_all_p, names_to="Method", values_to="Unpruned P-Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, 
         "Log-Unpruned" = -log10(`Unpruned P-Value`),
         Data  = rep(c("External", "Internal", "External", "Internal + External", "Internal",
                       "Internal + External", "Internal", "External", "Internal + External"), 100))

t1e_conf = pivot_longer(t1e_conf, prox_p:skat_all_p, names_to="Method", values_to="Pruned P-Value") %>%
  mutate("Log-Pruned" = -log10(`Pruned P-Value`))

results = cbind(t1e, t1e_conf[, 2:3])

results$Method = factor(results$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p",
                                                 "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p"),
                        labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
                                 "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O"))

results$Calculation = factor(results$Calculation)
results$Scenario = factor(results$Scenario)
results$MAF = factor(results$MAF)
results$Data = factor(results$Data, levels=c("Internal", "External", "Internal + External"))

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors = c("#009E73", "#0072B2", "#D55E00", "#CC79A7")


# P-value Distribution
ggplot(results,
       aes(x=`Log-Pruned`, y=`Log-Unpruned`, color = Method, shape = Method)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5) +
  geom_hline(yintercept=-log10(0.05), linetype=2, linewidth=0.5) +
  geom_vline(xintercept=-log10(0.05), linetype=2, linewidth=0.5) +
  scale_color_manual(values=colors) +
  facet_grid(~Data, space = "free_x", scales = "free_x") +
  labs(title = "Scatter Plot of P-values: 100% vs 99% Pruned", x = "100% vs 99% -log10pvalue", y = "100% vs 100% -log10pvalue") +
  theme_bw(base_size = 20)

# External only
# P-value Distribution
colors = c("#009E73", "#0072B2", "#CC79A7")

ggplot(subset(results, Data %in% "External"),
       aes(x=`Log-Pruned`, y=`Log-Unpruned`, color = Method, shape = Method)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.5) +
  geom_hline(yintercept=-log10(0.05), linetype=2, linewidth=0.5) +
  geom_vline(xintercept=-log10(0.05), linetype=2, linewidth=0.5) +
  scale_color_manual(values=colors) +
  scale_shape_manual(values=c(16, 4, 3))+
  facet_grid(~Data, space = "free_x", scales = "free_x") +
  #facet_grid(Data~factor(Method, levels = c("iECAT-O", "ProxECAT", "LogProx", "SKAT-O")), space = "free_x", scales = "free_x") +
  labs(title = "Scatter Plot of P-values: 99% vs 80% Pruned", x = "100% vs 80% -log10pvalue", y = "100% vs 99% -log10pvalue") +
  theme_bw(base_size = 20)


### Checks
results2 = subset(results, Data %in% "External")
results2 = results2 %>% filter(`Log-Unpruned` > -log10(0.05) & `Log-Pruned` < -log10(0.05))
