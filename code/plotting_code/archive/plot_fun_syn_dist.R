# Plot the distribution of fun and syn alleles for proxECAT
# Plot distribution of fun and syn variants and ratio of fun:syn variants
# for 100% AFR, 100% NFE, and admixed populations

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

dir = 'C:/Users/sagee/Documents/GitHub/masters_project/'
dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/'

# Read in csv file
# counts = read.table(paste0(dir, 'proxECAT_counts_expanded_', int_prune, "_v_", ext_prune, '.csv'), header = T, sep = ',')
counts = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_', int_prune, "_v_", ext_prune, '.csv'), header = T, sep = ',')
colnames(counts) <- c('Case-Fun (O)', 'Case-Syn (O)', 'Control-Fun (O)', 'Control-Syn (O)', 
                           'Control-Fun (E)', 'Control-Syn (E)', 'Control-Fun (O-E)',
                           'Control-Syn (O-E)', 'Ratio-Case', 'Ratio-Control', 'P-Value')

ctrls_fun_df <- rbind(data.frame(ctrl_fun = counts$`Control-Fun (O)`, Source = "Observed"),
                      data.frame(ctrl_fun = counts$`Control-Fun (E)`, Source = "Expected based on observed cases"))

ctrls_syn_df <- rbind(data.frame(ctrl_syn = counts$`Control-Syn (O)`, Source = "Observed"),
                      data.frame(ctrl_syn = counts$`Control-Syn (E)`, Source = "Expected based on observed cases"))

ctrls_df <- rbind(data.frame(ctrls = counts$`Control-Fun (O)`, Source = "Observed Functional Alleles"),
                  data.frame(ctrls = counts$`Control-Fun (E)`, Source = "Expected Functional Alleles based on observed cases"),
                  data.frame(ctrls = counts$`Control-Syn (O)`, Source = "Observed Synonymous Alleles"),
                  data.frame(ctrls = counts$`Control-Syn (E)`, Source = "Expected Synonymous Alleles based on observed cases"))

cases_df <- rbind(data.frame(cases = counts$`Case-Fun (O)`, Source = "Observed Functional Alleles in Cases"),
                  data.frame(cases = counts$`Case-Syn (O)`, Source = "Observed Synonymous Alleles in Cases"))

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors = c("#0072B2", "#CC79A7")
colors2 = c("#0072B2", "#BC9F4C", "#CC79A7", "#009E73")
colors_cases = c("#CC79A7", "#009E73")

# Make the density plots
p1 <- ggplot(ctrls_df, aes(x = ctrls, color = Source, fill = Source)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors2) +
        scale_fill_manual(values=colors2) +
        labs(title = paste0('Distribution of Functional and Synonymous Alleles (100% vs ', ext_prune, '% Pruned)\nPop: 100% NFE'), 
             x = "Number of Alleles", y = "Density") +
        theme_bw(base_size = 18)
p1
ggsave(file = paste0(dir_out, 'dist_fun_syn_', Pop2, '_100v', ext_prune, '.jpg'),
       plot = p1, height = 5, width = 12, units = 'in')

p2 <- ggplot(cases_df, aes(x = cases, color = Source, fill = Source)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors_cases) +
        scale_fill_manual(values=colors_cases) +
        labs(title = paste0('Distribution of Functional and Synonymous Alleles (100% Pruned Cases)\nPop: 100% NFE'), 
             x = "Number of Alleles", y = "Density") +
        theme_bw(base_size = 18)
p2
ggsave(file = paste0(dir_out, 'dist_fun_syn_', Pop2, '_100%_cases.jpg'),
       plot = p2, height = 5, width = 12, units = 'in')

# ggplot(ctrls_fun_df, aes(x = ctrl_fun, color = Source, fill = Source)) +
#   geom_density(alpha = 0.4) +
#   #scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#   scale_color_manual(values=colors) +
#   scale_fill_manual(values=colors) +
#   #facet_grid(Data~factor(Method, levels = c("iECAT-O", "ProxECAT", "LogProx", "SKAT-O")), space = "free_x", scales = "free_x") +
#   labs(title = "Distribution of Functional Alleles (100% vs 99% Pruned)", x = "Number of Functional Alleles", y = "Density") +
#   theme_bw(base_size = 20)
# 
# ggplot(ctrls_syn_df, aes(x = ctrl_syn, color = Source, fill = Source)) +
#   geom_density(alpha = 0.4) +
#   #scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
#   scale_color_manual(values=colors) +
#   scale_fill_manual(values=colors) +
#   #facet_grid(Data~factor(Method, levels = c("iECAT-O", "ProxECAT", "LogProx", "SKAT-O")), space = "free_x", scales = "free_x") +
#   labs(title = "Distribution of Synonymous Alleles (100% vs 99% Pruned)", x = "Number of Synonymous Alleles", y = "Density") +
#   theme_bw(base_size = 20)

################################################################################
# Read in files to compare dist of fun syn variants across populations
dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/gene_info/'
dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/Results/densityPlots/'

afr_100v100 = read.table(paste0(dir, 'fun_syn_counts_20K_AFR_100v100.txt'), header = T)
nfe_100v100 = read.table(paste0(dir, 'fun_syn_counts_20K_NFE_100v100.txt'), header = T)
afr_160v100v80 = read.table(paste0(dir, 'fun_syn_counts_20K_AFR_160v100v80.txt'), header = T)
nfe_160v100v80 = read.table(paste0(dir, 'fun_syn_counts_20K_NFE_160v100v80.txt'), header = T)
afr_nfe_160v100v80 = read.table(paste0(dir, 'fun_syn_counts_23K_AFR_NFE_160v100v80.txt'), header = T)

afr_100v100 = pivot_longer(afr_100v100, fun:ratio, names_to="variants", values_to="Count") %>%
  mutate(Pop = '20K_AFR', Pruning = '100v100', Source_Pruning = paste(variants, Pruning, sep = '_'), 
         Source_Pop = paste(variants, Pop, sep = '_'))
nfe_100v100 = pivot_longer(nfe_100v100, fun:ratio, names_to="variants", values_to="Count") %>%
  mutate(Pop = '20K_NFE', Pruning = '100v100', Source_Pruning = paste(variants, Pruning, sep = '_'), 
         Source_Pop = paste(variants, Pop, sep = '_'))
afr_160v100v80 = pivot_longer(afr_160v100v80, fun:ratio, names_to="variants", values_to="Count") %>%
  mutate(Pop = '20K_AFR', Pruning = '160v100v80', Source_Pruning = paste(variants, Pruning, sep = '_'), 
         Source_Pop = paste(variants, Pop, sep = '_'))
nfe_160v100v80 = pivot_longer(nfe_160v100v80, fun:ratio, names_to="variants", values_to="Count") %>%
  mutate(Pop = '20K_NFE', Pruning = '160v100v80', Source_Pruning = paste(variants, Pruning, sep = '_'), 
         Source_Pop = paste(variants, Pop, sep = '_'))
afr_nfe_160v100v80 = pivot_longer(afr_nfe_160v100v80, fun:ratio, names_to="variants", values_to="Count") %>%
  mutate(Pop = '23K_AFR_NFE', Pruning = '160v100v80', Source_Pruning = paste(variants, Pruning, sep = '_'), 
         Source_Pop = paste(variants, Pop, sep = '_'))

counts = rbind(afr_100v100, nfe_100v100, afr_160v100v80, nfe_160v100v80, afr_nfe_160v100v80)

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors = c("#0072B2", "#CC79A7")
colors2 = c("#0072B2", "#CC79A7", "#BC9F4C", "#009E73")
colors_cases = c("#CC79A7", "#009E73")
colors3 = c("#0072B2", "#BC9F4C", "#CC79A7", "#009E73", "#E69F00", "#999999")
colors_ratio = c("#0072B2", "#CC79A7", "#BC9F4C")

### 20K AFR 100v100 vs 160v100v80
p1 <- ggplot(counts %>% filter(Pop == '20K_AFR', variants != 'ratio'), aes(x = Count, color = Source, fill = Source)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors2) +
        scale_fill_manual(values=colors2) +
        labs(title = paste0('Distribution of Functional and Synonymous Variants by Pruning Pipeline (100% Pruned)\nPop: 100% AFR (20K)'), 
             x = "Number of Variants", y = "Density") +
        facet_wrap(~gene) +
        theme_bw(base_size = 18)
p1
ggsave(file = paste0(dir_out, 'dist_fun_syn_variants_20K_AFR_100v100_v_160v100v80.jpg'),
       plot = p1, height = 8, width = 15, units = 'in')

### 20K NFE 100v100 vs 160v100v80
p2 <- ggplot(counts %>% filter(Pop == '20K_NFE', variants != 'ratio'), aes(x = Count, color = Source_Pruning, fill = Source_Pruning)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors2) +
        scale_fill_manual(values=colors2) +
        labs(title = paste0('Distribution of Functional and Synonymous Variants by Pruning Pipeline (100% Pruned)\nPop: 100% NFE (20K)'), 
             x = "Number of Variants", y = "Density") +
        facet_wrap(~gene) +
        theme_bw(base_size = 18)
p2
ggsave(file = paste0(dir_out, 'dist_fun_syn_variants_20K_NFE_100v100_v_160v100v80.jpg'),
       plot = p2, height = 8, width = 15, units = 'in')

### 20K AFR vs NFE 160v100v80
p3 <- ggplot(counts %>% filter(Pop != '23K_AFR_NFE', variants != 'ratio', Pruning == '160v100v80'), aes(x = Count, color = Source_Pop, fill = Source_Pop)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors2) +
        scale_fill_manual(values=colors2) +
        labs(title = paste0('Distribution of Functional and Synonymous Variants by Population (100% Pruned)\nPipeline: 160v100v80'), 
             x = "Number of Variants", y = "Density") +
        facet_wrap(~gene) +
        theme_bw(base_size = 18)
p3
ggsave(file = paste0(dir_out, 'dist_fun_syn_variants_20K_AFR_v_NFE_160v100v80.jpg'),
       plot = p3, height = 8, width = 15, units = 'in')

### 20K AFR vs NFE 100v100
p4 <- ggplot(counts %>% filter(Pop != '23K_AFR_NFE', variants != 'ratio', Pruning == '100v100'), aes(x = Count, color = Source_Pop, fill = Source_Pop)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors2) +
        scale_fill_manual(values=colors2) +
        labs(title = paste0('Distribution of Functional and Synonymous Variants by Population (100% Pruned)\nPipeline: 100v100'), 
             x = "Number of Variants", y = "Density") +
        facet_wrap(~gene) +
        theme_bw(base_size = 18)
p4
ggsave(file = paste0(dir_out, 'dist_fun_syn_variants_20K_AFR_v_NFE_100v100.jpg'),
       plot = p4, height = 8, width = 15, units = 'in')

### 20K AFR vs 20K NFE vs 23K admixed ratio fun:syn 160v100v80
p5 <- ggplot(counts %>% filter(variants == 'ratio', Pruning == '160v100v80'), aes(x = Count, color = Source_Pop, fill = Source_Pop)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors_ratio) +
        scale_fill_manual(values=colors_ratio) +
        labs(title = paste0('Distribution of Ratio of Functional to Synonymous Variants by Population (100% Pruned)\nPipeline: 160v100v80'), 
             x = "Ratio", y = "Density") +
        facet_wrap(~gene) +
        theme_bw(base_size = 18)
p5
ggsave(file = paste0(dir_out, 'dist_ratio_fun_syn_variants_20K_AFR_v_20K_NFE_v_23K_AFR_NFE_160v100v80.jpg'),
       plot = p5, height = 8, width = 15, units = 'in')


### 20K AFR vs 20K NFE vs 23K admixed dist fun & syn 160v100v80
p6 <- ggplot(counts %>% filter(variants != 'ratio', Pruning == '160v100v80'), aes(x = Count, color = Source_Pop, fill = Source_Pop)) +
        geom_density(alpha = 0.4) +
        scale_color_manual(values=colors3) +
        scale_fill_manual(values=colors3) +
        labs(title = paste0('Distribution of Functional and Synonymous Variants by Population (100% Pruned)\nPipeline: 160v100v80'), 
             x = "Number of Variants", y = "Density") +
        facet_wrap(~gene) +
        theme_bw(base_size = 18)
p6
ggsave(file = paste0(dir_out, 'dist_fun_syn_variants_20K_AFR_v_20K_NFE_v_23K_AFR_NFE_160v100v80.jpg'),
       plot = p6, height = 8, width = 15, units = 'in')
