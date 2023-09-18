#####
# Summarize JUST type I error results for a particular iteration of the code
#####

# define function for calculating power
my.power <- function(values, alpha=0.05){
  values2 = values[!is.na(values)]
  sig = which(as.numeric(values2) <= alpha)
  out = length(sig)/length(values2)
  return(out)
}

Pop1 = "AFR"
Pop2 = "NFE"
scen = "s1"
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc10k'  #Number of common controls: 'cc5k' or 'cc10k'
int_prune = 100
ext_prune = 80

dir = 'C:/Users/sagee/Documents/HendricksLab/mastersProject/Results/cc10k/'

### Type 1 error

t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e_adj = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
t1e_homo = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)

t1e_out = apply(t1e, 2, my.power)
t1e_adj_out = apply(t1e_adj, 2, my.power)
t1e_homo_out = apply(t1e_homo, 2, my.power)

write.csv(t(as.data.frame(t1e_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(t(as.data.frame(t1e_adj_out)), paste0(dir, "T1e_all_adj_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(t(as.data.frame(t1e_homo_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)


# Recalculate p-values
dir = 'C:/Users/sagee/Documents/GitHub/masters_project/'
counts99 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_99.csv'), header = T, sep = ',')
counts95 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_95.csv'), header = T, sep = ',')
counts90 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_90.csv'), header = T, sep = ',')
counts80 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_80.csv'), header = T, sep = ',')
pow <- c()
pow <- c(pow, my.power(counts99$P.Value))
pow <- c(pow, my.power(counts95$P.Value))
pow <- c(pow, my.power(counts90$P.Value))
pow <- c(pow, my.power(counts80$P.Value))
