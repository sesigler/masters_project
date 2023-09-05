# Plot the distribution of fun and syn alleles for proxECAT

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

dir = 'C:/Users/sagee/OneDrive/Documents/GitHub/masters_project/'

# Read in csv file
counts = read.table(paste0(dir, 'proxECAT_counts_expanded_', int_prune, "_v_", ext_prune, '.csv'), header = T, sep = ',')
colnames(counts) <- c('Case-Fun (O)', 'Case-Syn (O)', 'Control-Fun (O)', 'Control-Syn (O)', 
                           'Control-Fun (E)', 'Control-Syn (E)', 'Control-Fun (O-E)',
                           'Control-Syn (O-E)', 'Ratio-Case', 'Ratio-Control', 'P-Value')

ctrls_fun_df <- rbind(data.frame(ctrl_fun = counts$`Control-Fun (O)`, type = "Observed"),
                      data.frame(ctrl_fun = counts$`Control-Fun (E)`, type = "Expected"))

ctrls_syn_df <- rbind(data.frame(ctrl_syn = counts$`Control-Syn (O)`, type = "Observed"),
                      data.frame(ctrl_syn = counts$`Control-Syn (E)`, type = "Expected"))

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors = c("#0072B2", "#CC79A7")

ggplot(ctrls_fun_df, aes(x = ctrl_fun, color = type, fill = type)) +
  geom_histogram(aes(y=after_stat(density)), position = "identity", alpha = 0.4, bins = 20) +
  geom_density(alpha=0.4) +
  #scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  #facet_grid(Data~factor(Method, levels = c("iECAT-O", "ProxECAT", "LogProx", "SKAT-O")), space = "free_x", scales = "free_x") +
  labs(title = "Distribution of Function Alleles (100% vs 99% Pruned)", x = "Number of Functional Alleles", y = "Density") +
  theme_bw(base_size = 20)

ggplot(ctrls_syn_df, aes(x = ctrl_syn, color = type, fill = type)) +
  geom_histogram(aes(y=after_stat(density)), position = "identity", alpha = 0.4, bins = 20) +
  geom_density(alpha=0.4) +
  #scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
  scale_color_manual(values=colors) +
  scale_fill_manual(values=colors) +
  #facet_grid(Data~factor(Method, levels = c("iECAT-O", "ProxECAT", "LogProx", "SKAT-O")), space = "free_x", scales = "free_x") +
  labs(title = "Distribution of Synonymous Alleles (100% vs 99% Pruned)", x = "Number of Synonymous Alleles", y = "Density") +
  theme_bw(base_size = 20)


