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
dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/AFR_NFE_pops/'
dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/Results/maf_mac_plots/'

# prop_ests = read.table(paste0(dir, "T1e_cc_prop_ests_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)
# counts = read.table(paste0(dir, "T1e_MACs_MAFs_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)
# Counts filtered to rare variants
counts_adj = read.table(paste0(dir, "T1e_MACs_MAFs_adj_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)
counts_unadj = read.table(paste0(dir, "T1e_MACs_MAFs_unadj_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)
counts_adj_unrounded = read.table(paste0(dir, "T1e_MACs_MAFs_adj_unrounded_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)

counts_unadj = counts_unadj %>% mutate(maf_case_unadj_cc = case_maf - unadj_cc_maf)
counts_adj = counts_adj %>% mutate(maf_case_adj_cc = case_maf - adj_cc_maf)
counts_adj_unrounded = counts_adj_unrounded %>% mutate(maf_case_adj_cc_unrounded = case_maf - adj_cc_maf_unrounded)

colnames(counts_adj) = colnames(counts_unadj) = colnames(counts_adj_unrounded) = c("case_mac", "case_maf", "cc_mac", "cc_maf", "gene", "id", "case_cc_maf_diff")
counts = rbind(counts_unadj, counts_adj, counts_adj_unrounded)
counts = counts %>% mutate(group = rep(c("Unadjusted", "Adjusted-Round then Sum", "Adjusted-Unrounded"), times=c(209282, 208571, 208571)))
counts$group = factor(counts$group, levels=c("Unadjusted", "Adjusted-Round then Sum", "Adjusted-Unrounded"))

same_id = intersect(counts_adj$id, counts_unadj$id)
same_id_2 = intersect(counts_unadj$id, counts_adj$id)

same_adj = counts_adj %>% filter(id %in% same_id)
same_unadj = counts_unadj %>% filter(id %in% same_id_2)

pairs(counts[, c(1, 3, 5)])
plot(prop_ests$maf_nfe, prop_ests$maf_afr)

# Boxplot of difference in MAF between cases and external controls
unadj = data.frame(group = 'Unadjusted External Controls', maf = counts_unadj$maf_case_unadj_cc)
adj = data.frame(group = 'Adjusted External Controls', maf = counts_adj$maf_case_adj_cc)
unrounded = data.frame(group = 'Adjusted Unrounded External Controls', maf = counts_adj_unrounded$maf_case_adj_cc_unrounded)
boxplot_data = rbind(unadj, adj, unrounded)
boxplot_data$group <- as.factor(boxplot_data$group)

p0 <- ggplot(boxplot_data, aes(x=group, y=maf, fill=group)) +
        geom_boxplot() +
        # scale_y_continuous(limits = c(-0.0001, 0.0001)) +
        labs(y = 'Case AF - External Control AF', x='Group', title = 'Scenario 2: Boxplots of Difference in MAFs between Cases and External Controls \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p0
ggsave(file = paste0(dir_out, 'boxplot_rare_case_cc_maf_diff_rounding_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
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
p3 <- ggplot(counts_adj, aes(x=adj_cc_mac, y=case_mac, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'Case MAC', x= 'External Control MAC (Adjusted)', title = 'Scenario 2: Scatter Plot of Case vs Adjusted External Control MACs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p3
ggsave(file = paste0(dir_out, 'rare_mac_case_v_adj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p3, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs adj cc MAFs
p4 <- ggplot(counts_adj, aes(x=adj_cc_maf, y=case_maf, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'Case AF', x= 'External Control AF (Adjusted)', title = 'Scenario 2: Scatter Plot of Case vs Adjusted External Control MAFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        # theme_bw(base_size = 18)
p4
ggsave(file = paste0(dir_out, 'rare_maf_case_v_adj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p4, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs unadj cc MACs
p5 <- ggplot(counts_unadj, aes(x=unadj_cc_mac, y=case_mac, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'Case MAC', x= 'External Control MAC (Unadjusted)', title = 'Scenario 2: Scatter Plot of Case vs Unadjusted External Control MACs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p5
ggsave(file = paste0(dir_out, 'rare_mac_case_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p5, height = 8, width = 16, units = 'in')

# Scatter plot by gene of case vs unadj cc MAFs
p6 <- ggplot(counts_unadj, aes(x=unadj_cc_maf, y=case_maf, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'Case AF', x= 'External Control AF (Unadjusted)', title = 'Scenario 2: Scatter Plot of Case vs Unadjusted External Control MAFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        # theme_bw(base_size = 18)
p6
ggsave(file = paste0(dir_out, 'rare_maf_case_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p6, height = 8, width = 16, units = 'in')

# Scatter plot by gene of adj vs unadj cc MACs
p7 <- ggplot(counts, aes(x=unadj_cc_mac, y=adj_cc_mac, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        # facet_wrap(~gene) +
        labs(y = 'External Control MAC (Adjusted)', x= 'External Control MAC (Unadjusted)', title = 'Scenario 2: Scatter Plot of Adjusted External Control vs Unadjusted External Control MACs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p7
ggsave(file = paste0(dir_out, 'mac_adj_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p7, height = 8, width = 16, units = 'in')

# Scatter plot by gene of adj vs unadj cc MAFs
p8 <- ggplot(counts, aes(x=unadj_cc_maf, y=adj_cc_maf, color=gene)) +
        geom_point() +
        geom_abline(slope = 1, intercept= 0) +
        facet_wrap(~gene) +
        labs(y = 'External Control AF (Adjusted)', x= 'External Control AF (Unadjusted)', title = 'Scenario 2: Scatter Plot of Adjusted External Control vs Unadjusted External Control MAFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme_bw(base_size = 18)
p8
ggsave(file = paste0(dir_out, 'maf_adj_v_unadj_cc_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p8, height = 8, width = 16, units = 'in')

cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette2 = c("#999999", "#BC9F4C", "#56B4E9", "#009E73", "#0072B2")
colors3 = c("#CC79A7", "#F0E442", "#56B4E9")
colors2 = c("#009E73", "#CC79A7")

# Scatter plot by gene of adj vs unadj cc MAFs
p9 <- ggplot(counts, aes(x=cc_maf, y=case_maf, color=group)) +
        geom_point(alpha=0.5) +
        geom_abline(slope = 1, intercept= 0) +
        scale_color_manual(values=colors3) +      
        facet_wrap(~gene) +
        labs(y = 'Case AF ', x= 'External Control AF', title = 'Scenario 2: Scatter Plot of Case vs External Control MAFs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
        theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        # theme_bw(base_size = 18)
p9
ggsave(file = paste0(dir_out, 'maf_case_v_cc_rounding_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p9, height = 8, width = 16, units = 'in')

# Scatter plot by gene of adj vs unadj cc MAFs
p10 <- ggplot(counts, aes(x=cc_mac, y=case_mac, color=group)) +
          geom_point(alpha=0.5) +
          geom_abline(slope = 1, intercept= 0) +
          scale_color_manual(values=colors3) +      
          facet_wrap(~gene) +
          labs(y = 'Case AC ', x= 'External Control AC', title = 'Scenario 2: Scatter Plot of Case vs External Control MACs \nCases: 100% AFR, External Controls: 80% AFR 20% NFE \nMAF=0.001') +
          theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
          # theme_bw(base_size = 18)
p10
ggsave(file = paste0(dir_out, 'mac_case_v_cc_rounding_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p10, height = 8, width = 16, units = 'in')

refs = read.table(paste0(dir, "ref_ac_af_", int_prune, "_v_", ext_prune, "_", Pop1, '_', Pop2, "_", scen, "_maf", maf, ".txt"), header = T)
refs$gene = factor(refs$gene, levels=c("ADGRE2", "ADGRE3", "ADGRE5", "CLEC17A", "DDX39A", "DNAJB1", 
                                       "GIPC1", "NDUFB7", "PKN1", "PTGER1", "TECR", "ZNF333"))
refs$fun = factor(refs$fun)
refs$common = factor(refs$common)

# Scatter plot ref AFs AFR vs NFE
p11 <- ggplot(refs, aes(x=maf_afr, y=maf_nfe, color=fun)) +
          geom_point(alpha=0.5) +
          geom_abline(slope = 1, intercept = 0) +
          scale_color_manual(values=colors2) +      
          # facet_wrap(~gene) +
          labs(y = 'NFE AF ', x= 'AFR AF', title = 'Scenario 2: Scatter Plot of AFR vs NFE AFs from Reference Data \nNref=500, Sim Reps=100 \nMAF=0.001') +
          theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
          # theme_bw(base_size = 18)
p11
ggsave(file = paste0(dir_out, 'ref_afr_vs_nfe_af_all_fun_syn_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p11, height = 8, width = 16, units = 'in')

# Scatter plot of 20K AFR vs 20K NFE AFs for the 100v100 common controls from 160v100v80
dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/checks/'
cc_afr = read.table(paste0(dir, "20K_", Pop1, "_cc_mac_maf_info_", int_prune, "_v_", ext_prune, "_maf", maf, ".txt"), header = T)
cc_nfe = read.table(paste0(dir, "20K_", Pop2, "_cc_mac_maf_info_", int_prune, "_v_", ext_prune, "_maf", maf, ".txt"), header = T)

merge <- merge(cc_afr, cc_nfe, by = "id", all = TRUE)

merge2 = merge %>% filter(maf.x <= maf & maf.y <= maf)

afr_data = cbind(merge2$mac.x, merge2$maf.x, )


p12 <- ggplot(merge2, aes(x=maf.x, y=maf.y)) +
          geom_point(alpha=0.5) +
          geom_abline(slope = 1, intercept = 0) +
          # scale_color_manual(values=colors2) +      
          # facet_wrap(~gene) +
          labs(y = 'NFE AF ', x= 'AFR AF', title = 'Scatter Plot of AFR vs NFE AFs from 20K Pop Data \nSim Reps=100 \nMAF=0.001') +
          theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
          # theme_bw(base_size = 18)
p12
ggsave(file = paste0(dir_out, 'ref_afr_vs_nfe_af_all_fun_syn_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p12, height = 8, width = 16, units = 'in')

################################################################################ 
#Plot case, adj and undj cc vs ref data
dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/checks/'
dir_out = 'C:/Users/sagee/Documents/GitHub/masters_project/Results/maf_mac_plots/'
cc_afr = read.table(paste0(dir, "ac_af_data_", int_prune, "_v_", ext_prune, "_", Pop1, "_", Pop2, "_", scen, "_maf", maf, ".txt"), header = T)

### Case vs Ref
# Case AF vs ref AFR AC
p13 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=case_maf)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Case AF ', x= 'Reference (AFR) AC', title = 'Scenario 2: Scatter Plot of Rare Case AF vs Reference AFR AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p13
ggsave(file = paste0(dir_out, 'rare_case_AF_vs_refAFR_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p13, height = 8, width = 16, units = 'in')

# Case AF vs ref NFE AC
p14 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=case_maf)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Case AF ', x= 'Reference (NFE) AC', title = 'Scenario 2: Scatter Plot of Rare Case AF vs Reference NFE AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p14
ggsave(file = paste0(dir_out, 'rare_case_AF_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p14, height = 8, width = 16, units = 'in')

# Case AC vs ref AFR AC
p15 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=case_mac)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Case AC ', x= 'Reference (AFR) AC', title = 'Scenario 2: Scatter Plot of Rare Case AC vs Reference AFR AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p15
ggsave(file = paste0(dir_out, 'rare_case_AC_vs_refAFR_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p15, height = 8, width = 16, units = 'in')

# Case AC vs ref NFE AC
p16 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=case_mac)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Case AC ', x= 'Reference (NFE) AC', title = 'Scenario 2: Scatter Plot of Rare Case AC vs Reference NFE AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p16
ggsave(file = paste0(dir_out, 'rare_case_AC_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p16, height = 8, width = 16, units = 'in')

### Adj CC vs Ref
# adj CC AF vs ref AFR AC
p17 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=adj_maf)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Adjusted External Control AF', x= 'Reference (AFR) AC', title = 'Scenario 2: Scatter Plot of Rare Adjusted External Control AF vs Reference AFR AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p17
ggsave(file = paste0(dir_out, 'rare_adj_cc_AF_vs_refAFR_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p17, height = 8, width = 16, units = 'in')

# adj CC AF vs ref NFE AC
p18 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=adj_maf)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Adjusted External Control AF', x= 'Reference (NFE) AC', title = 'Scenario 2: Scatter Plot of Rare Adjusted External Control AF vs Reference NFE AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p18
ggsave(file = paste0(dir_out, 'rare_adj_cc_AF_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p18, height = 8, width = 16, units = 'in')

# adj CC AC vs ref AFR AC
p19 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=adj_mac)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Adjusted External Control AC', x= 'Reference (AFR) AC', title = 'Scenario 2: Scatter Plot of Rare Adjusted External Control AC vs Reference AFR AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p19
ggsave(file = paste0(dir_out, 'rare_adj_cc_AC_vs_refAFR_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p19, height = 8, width = 16, units = 'in')

# adj CC AC vs ref NFE AC
p20 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=adj_mac)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Adjusted External Control AC', x= 'Reference (NFE) AC', title = 'Scenario 2: Scatter Plot of Rare Adjusted External Control AC vs Reference NFE AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p20
ggsave(file = paste0(dir_out, 'rare_adj_cc_AC_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p20, height = 8, width = 16, units = 'in')

### Unadj CC vs Ref
# unadj CC AF vs ref AFR AC
p21 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=maf)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Unadjusted External Control AF', x= 'Reference (AFR) AC', title = 'Scenario 2: Scatter Plot of Rare Unadjusted External Control AF vs Reference AFR AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p21
ggsave(file = paste0(dir_out, 'rare_unadj_cc_AF_vs_refAFR_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p21, height = 8, width = 16, units = 'in')

# unadj CC AF vs ref NFE AC
p22 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=maf)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Unadjusted External Control AF', x= 'Reference (NFE) AC', title = 'Scenario 2: Scatter Plot of Rare Unadjusted External Control AF vs Reference NFE AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p22
ggsave(file = paste0(dir_out, 'rare_unadj_cc_AF_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p22, height = 8, width = 16, units = 'in')

# unadj CC AC vs ref AFR AC
p23 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=mac)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Unadjusted External Control AC', x= 'Reference (AFR) AC', title = 'Scenario 2: Scatter Plot of Rare Unadjusted External Control AC vs Reference AFR AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p23
ggsave(file = paste0(dir_out, 'rare_unadj_cc_AC_vs_refAFR_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p23, height = 8, width = 16, units = 'in')

# unadj CC AC vs ref NFE AC
p24 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=mac)) +
        geom_point(alpha=0.5) +
        # geom_abline(slope = 1, intercept = 0) +
        # scale_color_manual(values=colors2) +      
        facet_wrap(~gene) +
        labs(y = 'Unadjusted External Control AC', x= 'Reference (NFE) AC', title = 'Scenario 2: Scatter Plot of Rare Unadjusted External Control AC vs Reference NFE AC from 160v100v80 Data \nSim Reps=100, Nref=500, Ncase=5000 \nMAF=0.001') +
        # theme(axis.text.x = element_text(angle = 35, hjust=0.65)) 
        theme_bw(base_size = 18)
p24
ggsave(file = paste0(dir_out, 'rare_unadj_cc_AC_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = p24, height = 8, width = 16, units = 'in')

### Merge plots into one plot
# AF vs refs
p13 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=case_maf)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Case AF ', x= 'Reference (AFR) AC')

p14 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=case_maf)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Case AF ', x= 'Reference (NFE) AC')

p17 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=adj_maf)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Adjusted External Control AF', x= 'Reference (AFR) AC')

p18 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=adj_maf)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Adjusted External Control AF', x= 'Reference (NFE) AC')

p21 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=maf)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Unadjusted External Control AF', x= 'Reference (AFR) AC')

p22 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=maf)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Unadjusted External Control AF', x= 'Reference (NFE) AC')

af_v_ref <- ggarrange(p13, p14, p17, p18, p21, p22, ncol=2, nrow=3)
af_v_ref
ggsave(file = paste0(dir_out, 'rare_AF_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = af_v_ref, height = 8, width = 16, units = 'in')


# AC vs refs
p15 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=case_mac)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Case AC ', x= 'Reference (AFR) AC') 

p16 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=case_mac)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Case AC ', x= 'Reference (NFE) AC')

p19 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=adj_mac)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Adjusted External Control AC', x= 'Reference (AFR) AC')

p20 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=adj_mac)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Adjusted External Control AC', x= 'Reference (NFE) AC')

p23 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_afr, y=mac)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Unadjusted External Control AC', x= 'Reference (AFR) AC') 

p24 <- ggplot(cc_afr %>% filter(gene == "ZNF333"), aes(x=mac_nfe, y=mac)) +
  geom_point(alpha=0.5) +
  facet_wrap(~gene) +
  labs(y = 'Unadjusted External Control AC', x= 'Reference (NFE) AC')

ac_v_ref <- ggarrange(p15, p16, p19, p20, p23, p24, ncol=2, nrow=3)
ac_v_ref
ggsave(file = paste0(dir_out, 'rare_AC_vs_refNFE_AC_', Pop1, '_', Pop2, '_', scen, '_', folder, '_', int_prune, '_v_', ext_prune, '_maf', maf, '.jpg'), 
       plot = ac_v_ref, height = 8, width = 16, units = 'in')

# Histogram of Case AF in variants with MAC 0 in AFR ref
p_hist_caseAF_afr0 <- ggplot(cc_afr %>% filter(gene == "ZNF333" & mac_afr == 0), aes(x=case_maf)) +
                        geom_histogram(bins=20, color = "black", fill = "white") +
                        facet_wrap(~gene)
p_hist_caseAF_afr0

p_hist_adjccAF_afr0 <- ggplot(cc_afr %>% filter(gene == "ZNF333" & mac_afr == 0), aes(x=adj_maf)) +
                        geom_histogram(bins=20, color = "black", fill = "white") +
                        facet_wrap(~gene)
p_hist_adjccAF_afr0

p_hist_unadjccAF_afr0 <- ggplot(cc_afr %>% filter(gene == "ZNF333" & mac_afr == 0), aes(x=maf)) +
                          geom_histogram(bins=20, color = "black", fill = "white") +
                          facet_wrap(~gene)
p_hist_unadjccAF_afr0

p_hist_caseAF_afr1 <- ggplot(cc_afr %>% filter(gene == "ZNF333" & mac_afr == 1), aes(x=case_maf)) +
                        geom_histogram(bins=20, color = "black", fill = "white") +
                        facet_wrap(~gene)
p_hist_caseAF_afr1

p_hist_adjccAF_afr1 <- ggplot(cc_afr %>% filter(gene == "ZNF333" & mac_afr == 1), aes(x=adj_maf)) +
                        geom_histogram(bins=20, color = "black", fill = "white") +
                        facet_wrap(~gene)
p_hist_adjccAF_afr1

p_hist_unadjccAF_afr1 <- ggplot(cc_afr %>% filter(gene == "ZNF333" & mac_afr == 1), aes(x=maf)) +
                          geom_histogram(bins=20, color = "black", fill = "white") +
                          facet_wrap(~gene)
p_hist_unadjccAF_afr1
