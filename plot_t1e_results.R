### Used to plot type I error results ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

Pop1 = 'AFR'
Pop2 = 'NFE'
scen = 's1'
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc5k'  #Number of common controls: 'cc5k' or 'cc10k'
int_prune = 80
ext_prune = 80

dir = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/Results/cc5k/'

# read in the results
# Proportion Estimates 
# t1e_case_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_case_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_cc_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_cc_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_int_prop_ests = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_int_prop_ests_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

# t1e
t1e_all = read.csv(paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
t1e_all_adj = read.csv(paste0(dir, "T1e_all_adj_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
t1e_all_homo = read.csv(paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)

# puts in a format for ggplot
t1e_all = pivot_longer(t1e_all, prox_p:skat_all_p, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "Admixed")
t1e_all_adj = pivot_longer(t1e_all_adj, prox_p_adj:iecat_p_adj, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "Admixed")
t1e_all_homo = pivot_longer(t1e_all_homo, prox_p:skat_all_p, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = scen, MAF = maf, Pop = "Homogeneous")

results = rbind(t1e_all_homo, t1e_all, t1e_all_adj)

results2 = results %>% mutate(Data = c("External", "Internal", "External", "Internal + External", "Internal", 
                                       "Internal + External", "Internal", "External", "Internal + External",
                                       "External", "Internal", "External", "Internal + External", "Internal", 
                                       "Internal + External", "Internal", "External", "Internal + External",
                                       "External", "External", "Internal + External", "Internal + External")) %>%
  mutate(MACs = rep(c("Unadjusted", "Adjusted"), c(18, 4))) %>%
  mutate(Configuration = rep(c("Homogeneous-Unadjusted", "Admixed-Unadjusted", "Admixed-Adjusted"), c(9, 9, 4)))

results2$Method = factor(results2$Method, levels=c("prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
                                                   "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p",
                                                   "prox_p", "prox_int_p", "prox2_p", "prox2_all_p", "prox2_int_p", 
                                                   "iecat_p", "skat_int_p", "skat_ext_p", "skat_all_p",
                                                   "prox_p_adj", "prox2_p_adj", "prox2_all_p_adj", "iecat_p_adj"),
                         labels=c("ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
                                  "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O",
                                  "ProxECAT", "ProxECAT", "LogProx", "LogProx", "LogProx", 
                                  "iECAT-O", "SKAT-O", "SKAT-O", "SKAT-O",
                                  "ProxECAT", "LogProx", "LogProx", "iECAT-O"))


#results2$MAF = factor(results2$MAF, levels = c(0.01, 0.001))
results2$Calculation = factor(results2$Calculation)
results2$Scenario = factor(results2$Scenario)
results2$MAF = factor(results2$MAF)
results2$Pop = factor(results2$Pop, levels=c("Admixed", "Homogeneous"))
results2$Data = factor(results2$Data, levels=c("Internal", "External", "Internal + External"))
results2$MACs = factor(results2$MACs, levels=c("Unadjusted", "Adjusted"))
results2$Configuration = factor(results2$Configuration, levels=c("Homogeneous-Unadjusted", "Admixed-Unadjusted", "Admixed-Adjusted"))

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
colors = c("#009E73", "#0072B2", "#D55E00")

# Controls everything in the graph
# base_size = 25

# Type I error graph
ggplot(results2 %>% filter(Calculation=="Type I Error", MAF==0.001, Scenario==scen),  
       aes(x=Method, y=Value, color=Configuration)) +
  geom_point(aes(shape=Configuration), size=5, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
  scale_y_continuous(limits=c(0, 1)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
  scale_color_manual(values=colors) +
  facet_grid(~Data, scales="free", space="free") +
  labs(y='Type I Error', x='Method', title='Scenario 1: Type I Error 80% vs 80% (5k cc) \nMAF=0.001') +
  theme_bw(base_size = 23)


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
