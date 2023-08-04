### Used to plot Master's project results ###

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(binom)
library(data.table) # for fread

Pop1 = 'AFR'
Pop2 = 'NFE'
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc10k'  #Number of common controls: 'cc5k' or 'cc10k'

dir = 'C:/Users/sagee/OneDrive/Documents/HendricksLab/mastersProject/Results/cc10k/'

# read in the results
# Proportion Estimates 
t1e.case.prop.ests.s1 = read.table(paste0(dir, "T1e_Confounding_case_prop_ests_s1_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e.case.prop.ests.s2 = read.table(paste0(dir, "T1e_Confounding_case_prop_ests_s2_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e.cc.prop.ests.s1 = read.table(paste0(dir, "T1e_Confounding_cc_prop_ests_s1_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e.cc.prop.ests.s2 = read.table(paste0(dir, "T1e_Confounding_cc_prop_ests_s2_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e.int.prop.ests.s1 = read.table(paste0(dir, "T1e_Confounding_int_prop_ests_s1_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e.int.prop.ests.s2 = read.table(paste0(dir, "T1e_Confounding_int_prop_ests_s2_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

power.case.prop.ests.s1 = read.table(paste0(dir, "Power_Confounding_case_prop_ests_s1_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
power.case.prop.ests.s2 = read.table(paste0(dir, "Power_Confounding_case_prop_ests_s2_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
power.cc.prop.ests.s1 = read.table(paste0(dir, "Power_Confounding_cc_prop_ests_s1_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
power.cc.prop.ests.s2 = read.table(paste0(dir, "Power_Confounding_cc_prop_ests_s2_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
power.int.prop.ests.s1 = read.table(paste0(dir, "Power_Confounding_int_prop_ests_s1_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
power.int.prop.ests.s2 = read.table(paste0(dir, "Power_Confounding_int_prop_ests_s2_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)

# t1e
t1e.all.s1 = read.csv(paste0(dir, "T1e_Confounding_all_s1_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
t1e.all.s1.adj = read.csv(paste0(dir, "T1e_Confounding_all_s1_adj_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
t1e.all.s2 = read.csv(paste0(dir, "T1e_Confounding_all_s2_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
t1e.all.s2.adj = read.csv(paste0(dir, "T1e_Confounding_all_s2_adj_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), header=T)
# power
power.all.s1 = read.csv(paste0(dir, "Power_Confounding_all_s1_", Pop1, '-', Pop2, "_maf", maf, ".csv"), header=T)
power.all.s1.adj = read.csv(paste0(dir, "Power_Confounding_all_s1_adj_", Pop1, '-', Pop2, "_maf", maf, ".csv"), header=T)
power.all.s2 = read.csv(paste0(dir, "Power_Confounding_all_s2_", Pop1, '-', Pop2, "_maf", maf, ".csv"), header=T)
power.all.s2.adj = read.csv(paste0(dir, "Power_Confounding_all_s2_adj_", Pop1, '-', Pop2, "_maf", maf, ".csv"), header=T)

# puts in a format for ggplot
# t1e
# type1.conf.nfe = pivot_longer(type1.conf.nfe, prox.p:skat.all2.p, names_to="Method", values_to="Power") %>%
#   mutate(Scenario = "Type I Error", Confounding = "Yes", Population = "NFE", MAF = maf)
t1e.all.s1 = pivot_longer(t1e.all.s1, prox.p.s1:iecat.p.s1, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = "S1", MAF = maf)
t1e.all.s1.adj = pivot_longer(t1e.all.s1.adj, prox.p.s1.adj:iecat.p.s1.adj, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = "S1", MAF = maf)
t1e.all.s2 = pivot_longer(t1e.all.s2, prox.p.s2:iecat.p.s2, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = "S2", MAF = maf)
t1e.all.s2.adj = pivot_longer(t1e.all.s2.adj, prox.p.s2.adj:iecat.p.s2.adj, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Type I Error", Scenario = "S2", MAF = maf)

# power
power.all.s1 = pivot_longer(power.all.s1, prox.p.s1:iecat.p.s1, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Power", Scenario = "S1", MAF = maf)
power.all.s1.adj = pivot_longer(power.all.s1.adj, prox.p.s1.adj:iecat.p.s1.adj, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Power", Scenario = "S1", MAF = maf)
power.all.s2 = pivot_longer(power.all.s2, prox.p.s2:iecat.p.s2, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Power", Scenario = "S2", MAF = maf)
power.all.s2.adj = pivot_longer(power.all.s2.adj, prox.p.s2.adj:iecat.p.s2.adj, names_to="Method", values_to="Value") %>%
  mutate(Calculation = "Power", Scenario = "S2", MAF = maf)

# results = rbind(t1e.all.s1, t1e.all.s1.adj, t1e.all.s2, t1e.all.s2.adj,
#                 power.all.s1, power.all.s1.adj, power.all.s2, power.all.s2.adj)
results = rbind(t1e.all.s1, t1e.all.s1.adj, t1e.all.s2, t1e.all.s2.adj)

# results2 = results %>% filter(Method %in% c("prox.p", "prox.int.p", "prox2.p", "prox2.int.p", "prox2.all.p", "iecat.p")) %>%
#   mutate(Data = rep(c("External", "External", "Internal + External", "Internal + External", "Internal", "External"), 8))

results2 = results %>% mutate(Data = rep(c("External", "External", "Internal + External", "Internal + External"), 4)) %>%
  mutate(MACs = rep(c("Unadjusted", "Unadjusted","Unadjusted","Unadjusted",
                            "Adjusted", "Adjusted", "Adjusted", "Adjusted"), 2))

results2$Method = factor(results2$Method, levels=c("prox.p.s1", "prox2.p.s1", "prox2.all.p.s1", "iecat.p.s1",
                                                   "prox.p.s1.adj", "prox2.p.s1.adj",  "prox2.all.p.s1.adj", "iecat.p.s1.adj",
                                                   "prox.p.s2", "prox2.p.s2", "prox2.all.p.s2", "iecat.p.s2",
                                                   "prox.p.s2.adj", "prox2.p.s2.adj", "prox2.all.p.s2.adj", "iecat.p.s2.adj"),
                         labels=c("ProxECAT", "LogProx", "LogProx", "iECAT-O",
                                  "ProxECAT", "LogProx", "LogProx", "iECAT-O",
                                  "ProxECAT", "LogProx", "LogProx", "iECAT-O",
                                  "ProxECAT", "LogProx", "LogProx", "iECAT-O"))

results2$Data = factor(results2$Data, levels=c("External", "Internal + External"))
#results2$MAF = factor(results2$MAF, levels = c(0.01, 0.001))
results2$MAF = factor(results2$MAF)
results2$Scenario = factor(results2$Scenario, levels = c("S1", "S2"))

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
colors = c("#0072B2", "#D55E00")

# Controls everything in the graph
# base_size = 25

# Type I error graph
# Scenario 1
ggplot(results2 %>% filter(Calculation=="Type I Error", MAF==0.001, Scenario=="S1"),  
       aes(x=Method, y=Value, color=MACs)) +
  geom_point(aes(shape=MACs), size=5, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
  scale_y_continuous(limits=c(0, 1)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
  scale_color_manual(values=colors) +
  facet_grid(~Data, scales="free", space="free") +
  labs(y='Type I Error', x='Method', title='Scenario 1: Type I Error (10k cc) \nMAF=0.001') +
  theme_bw(base_size = 23)

# Scenario 2
ggplot(results2 %>% filter(Calculation=="Type I Error", MAF==0.001, Scenario=="S2"),  
       aes(x=Method, y=Value, color=MACs)) +
  geom_point(aes(shape=MACs), size=5, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.05, linetype=2, linewidth=1.5) +
  scale_y_continuous(limits=c(0, 1)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), linewidth=1.5, width=.2, position=position_dodge(width=0.5)) +
  scale_color_manual(values=colors) +
  facet_grid(~Data, scales="free", space="free") +
  labs(y='Type I Error', x='Method', title='Scenario 2: Type I Error (10k cc) \nMAF=0.001') +
  theme_bw(base_size = 23)

# results2$Power[which(results2$Method=="SKAT-O" & results2$Scenario=="Power" & results2$Confounding=="Yes" & results2$Data=="External")] = NA
# results2$Lower[which(results2$Method=="SKAT-O" & results2$Scenario=="Power" & results2$Confounding=="Yes" & results2$Data=="External")] = NA
# results2$Upper[which(results2$Method=="SKAT-O" & results2$Scenario=="Power" & results2$Confounding=="Yes" & results2$Data=="External")] = NA

#%>% filter(!(Method=="SKAT-O" & Data=="External"))
# Power graph
# Scenario 1
ggplot(results2 %>% filter(Calculation=="Power", MAF==0.01, Scenario=="S1"),  
       aes(x=Method, y=Value, color=MACs)) +
  geom_point(aes(shape=MACs), size=5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.5, width=.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.8, linetype=2, size=1.5) +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values=colors) +
  facet_grid(~Data, scales="free", space="free") +
  labs(y='Power', x='Method', title='Scenario 1: Power \nMAF=0.01') +
  theme_bw(base_size=23)
  #theme(axis.text.x=element_text(angle=30, vjust=1, hjust=1))

# Scenario 2
ggplot(results2 %>% filter(Calculation=="Power", MAF==0.01, Scenario=="S2"),  
       aes(x=Method, y=Value, color=MACs)) +
  geom_point(aes(shape=MACs), size=5, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin=Lower, ymax=Upper), size=1.5, width=.2, position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0.8, linetype=2, size=1.5) +
  scale_y_continuous(limits=c(0,1)) +
  scale_color_manual(values=colors) +
  facet_grid(~Data, scales="free", space="free") +
  labs(y='Power', x='Method', title='Scenario 2: Power \nMAF=0.01') +
  theme_bw(base_size=23)

# Proportion estimate plots
# t1e Scenario 1
par(mfrow=c(2,3))
plot(t1e.cc.prop.ests.s1$maf_afr, t1e.case.prop.ests.s1$maf_afr, 
     xlab = "AFR AF External Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(t1e.int.prop.ests.s1$maf_afr, t1e.case.prop.ests.s1$maf_afr, 
     xlab = "AFR AF Internal Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(t1e.cc.prop.ests.s1$maf_afr, t1e.int.prop.ests.s1$maf_afr, 
     xlab = "AFR AF External Controls", ylab = "AFR AF Internal Controls")
abline(a=0, b=1)
mtext("Scenario 1 Dataset Proportion Estimates-Type I error", side = 3, line = -3, outer = TRUE)

# t1e Scenario 2
plot(t1e.cc.prop.ests.s2$maf_afr, t1e.case.prop.ests.s2$maf_afr, xlim = c(0.8, 1),
     ylim = c(0.8, 1), 
     xlab = "AFR AF External Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(t1e.int.prop.ests.s2$maf_afr, t1e.case.prop.ests.s2$maf_afr, 
     xlab = "AFR AF Internal Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(t1e.cc.prop.ests.s2$maf_afr, t1e.int.prop.ests.s2$maf_afr, xlim = c(0.8, 1),
     ylim = c(0.8, 1), 
     xlab = "AFR AF External Controls", ylab = "AFR AF Internal Controls")
abline(a=0, b=1)
mtext("Scenario 2 Dataset Proportion Estimates-Type I error", side = 3, line = -29, outer = TRUE)

# Power Scenario 1
par(mfrow=c(2,3))
plot(power.cc.prop.ests.s1$maf_afr, power.case.prop.ests.s1$maf_afr, 
     xlab = "AFR AF External Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(power.int.prop.ests.s1$maf_afr, power.case.prop.ests.s1$maf_afr, 
     xlab = "AFR AF Internal Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(power.cc.prop.ests.s1$maf_afr, power.int.prop.ests.s1$maf_afr, 
     xlab = "AFR AF External Controls", ylab = "AFR AF Internal Controls")
abline(a=0, b=1)
mtext("Scenario 1 Dataset Proportion Estimates-Power", side = 3, line = -3, outer = TRUE)

# Power Scenario 2
plot(power.cc.prop.ests.s2$maf_afr, power.case.prop.ests.s2$maf_afr, xlim = c(0.8, 1),
     ylim = c(0.8, 1), 
     xlab = "AFR AF External Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(power.int.prop.ests.s2$maf_afr, power.case.prop.ests.s2$maf_afr,
     xlab = "AFR AF Internal Controls", ylab = "AFR AF Cases")
abline(a=0, b=1)

plot(power.cc.prop.ests.s2$maf_afr, power.int.prop.ests.s2$maf_afr, xlim = c(0.8, 1),
     ylim = c(0.8, 1), 
     xlab = "AFR AF External Controls", ylab = "AFR AF Internal Controls")
abline(a=0, b=1)
mtext("Scenario 2 Dataset Proportion Estimates-Power", side = 3, line = -29, outer = TRUE)


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
