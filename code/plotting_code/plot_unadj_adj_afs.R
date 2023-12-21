# Plot boxplot of difference in AFs between cases and external controls
# Plot scatterplot of cases, unadjusted external controls and adjusted external controls
# for admixed populations, both scenarios

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread
library(GGally)

Pop1 = 'AFR'
Pop2 = 'NFE'
scen = 's2'
maf = 0.001
folder = '160v100v80'
int_prune = 100
ext_prune = 100

dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/AFR_NFE_pops/af_checks/'
dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/Results/maf_mac_plots/'

prop_ests = read.table(paste0(dir, "T1e_cc_prop_ests_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)
counts = read.table(paste0(dir, "T1e_MACs_MAFs_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)

counts = counts %>% mutate(maf_case_unadj_cc = case_maf - unadj_cc_maf,
                           maf_case_adj_cc = case_maf - adj_cc_maf)

pairs(counts[, c(1, 3, 5)])
plot(prop_ests$maf_nfe, prop_ests$maf_afr)

# Boxplot of difference in MAF between cases and external controls
unadj = data.frame(group = 'Unadjusted External Controls', value = counts$maf_case_unadj_cc)
adj = data.frame(group = 'Adjusted External Controls', value = counts$maf_case_adj_cc)
boxplot_data = rbind(unadj, adj)
boxplot_data$group <- as.factor(boxplot_data$group)

p0 <- ggplot(boxplot_data, aes(x=group, y=abs(value), fill=group)) +
  geom_boxplot() +
  labs(y = '|Case AF - External Control AF|', x='Group', title = 'Sceanrio 2: Boxplots of Difference in MAFs between Cases and External Controls \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
  theme_bw(base_size = 18)
p0
ggsave(file = paste0(dir_out, 'boxplot_case_cc_maf_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p0, height = 8, width = 16, units = 'in')

# Scatter plot matrix of AFs
p1 <- ggpairs(counts, columns = c(2, 4, 6), aes(color=gene),
              title = "Scenario 2: Scatter Plot Matrix of MAFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001")
p1
ggsave(file = paste0(dir_out, 'maf_matrix_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p1, height = 8, width = 16, units = 'in')

# Scatter plot matrix of ACs
p2 <- ggpairs(counts, columns = c(1, 3, 5), aes(color=gene),
              title = "Scenario 2: Scatter Plot Matrix of MACs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001")
p2
ggsave(file = paste0(dir_out, 'mac_matrix_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p2, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs adj cc MACs
p3 <- ggplot(counts, aes(x=adj_cc_mac, y=case_mac, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'Case MAC', x= 'External Control MAC (Adjusted)', title = 'Scenario 2: Scatter Plot of Case vs Adjusted External Control AFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p3
ggsave(file = paste0(dir_out, 'mac_case_v_adj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p3, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs adj cc MAFs
p4 <- ggplot(counts, aes(x=adj_cc_maf, y=case_maf, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        # facet_wrap(~gene) +
        labs(y = 'Case AF', x= 'External Control AF (Adjusted)', title = 'Scatter Plots of Case vs Adjusted External Control AFs by Gene \nScenario 2: Cases 100% AFR, External Controls 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p4
ggsave(file = paste0(dir_out, 'maf_case_v_adj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p4, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs unadj cc MACs
p5 <- ggplot(counts, aes(x=unadj_cc_mac, y=case_mac, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        # facet_wrap(~gene) +
        labs(y = 'Case MAC', x= 'External Control MAC (Unadjusted)', title = 'Scenario 2: Scatter Plot of Case vs Unadjusted External Control AFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p5
ggsave(file = paste0(dir_out, 'mac_case_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p5, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs unadj cc MAFs
p6 <- ggplot(counts, aes(x=unadj_cc_maf, y=case_maf, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        # facet_wrap(~gene) +
        labs(y = 'Case AF', x= 'External Control AF (Unadjusted)', title = 'Scatter Plots of Case vs Unadjusted External Control AFs by Gene \nScenario 2: Cases 100% AFR, External Controls 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p6
ggsave(file = paste0(dir_out, 'maf_case_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p6, height = 8, width = 16, units = 'in')

# Scatter plot by gene of adj vs unadj cc MACs
p7 <- ggplot(counts, aes(x=unadj_cc_mac, y=adj_cc_mac, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        # facet_wrap(~gene) +
        labs(y = 'External Control MAC (Adjusted)', x= 'External Control MAC (Unadjusted)', title = 'Scenario 2: Scatter Plot of Adjusted External Control vs Unadjusted External Control AFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p7
ggsave(file = paste0(dir_out, 'mac_adj_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p7, height = 8, width = 16, units = 'in')

# Scatter plot by gene of adj vs unadj cc MAFs
p8 <- ggplot(counts, aes(x=unadj_cc_maf, y=adj_cc_maf, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'External Control AF (Adjusted)', x= 'External Control AF (Unadjusted)', title = 'Scatter Plots of Adjusted External Control vs Unadjusted External Control AFs by Gene \nScenario 2: Cases 100% AFR, External Controls 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p8
ggsave(file = paste0(dir_out, 'maf_adj_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p8, height = 8, width = 16, units = 'in')




