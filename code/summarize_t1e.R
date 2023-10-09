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
int_prune = 80
ext_prune = 80
pruning = "pruneSepRaresim" #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '100v80'

dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', folder, '/')

### Type 1 error

# t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_adj = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_homo = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_skat = read.table(paste0(dir, "T1e_skat_syn_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_og_hap = read.table(paste0(dir, "T1e_OG_hap_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
t1e_pruning = read.table(paste0(dir, "T1e_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_fixed_proxecat_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)

# t1e_out = apply(t1e, 2, my.power)
# t1e_adj_out = apply(t1e_adj, 2, my.power)
# t1e_homo_out = apply(t1e_homo, 2, my.power)
# t1e_skat_out = apply(t1e_skat, 2, my.power)
# t1e_og_hap_out = apply(t1e_og_hap, 2, my.power)
t1e_pruning_out = apply(t1e_pruning, 2, my.power)

# write.csv(t(as.data.frame(t1e_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_adj_out)), paste0(dir, "T1e_all_adj_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_homo_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_skat_out)), paste0(dir, "T1e_all_skat_syn_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_og_hap_out)), paste0(dir, "T1e_all_OG_hap_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)



# Recalculate p-values proxECAT only data
dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/'
counts99 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_99.csv'), header = T, sep = ',')
counts95 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_95.csv'), header = T, sep = ',')
counts90 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_90.csv'), header = T, sep = ',')
counts80 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_80.csv'), header = T, sep = ',')
pow <- data.frame(matrix(nrow = 1, ncol = 4))
colnames(pow) <- c('t1e_99', 't1e_95', 't1e_90', 't1e_80')
pow[1, 1] <- my.power(counts99$P.Value)
pow[1, 2] <- my.power(counts95$P.Value)
pow[1, 3] <- my.power(counts90$P.Value)
pow[1, 4] <- my.power(counts80$P.Value)

write.csv(pow, paste0(dir, "T1e_NFE_99-80_maf", maf, ".csv"), quote=F, row.names=F)

counts100v100 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_100.csv'), header = T, sep = ',')
counts99v99 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_99_v_99.csv'), header = T, sep = ',')
counts95v95 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_95_v_95.csv'), header = T, sep = ',')
counts90v90 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_90_v_90.csv'), header = T, sep = ',')
counts80v80 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_80_v_80.csv'), header = T, sep = ',')

t1e_ext_v_ext <- data.frame(matrix(nrow = 1, ncol = 5))
colnames(t1e_ext_v_ext) <- c('t1e_100v100', 't1e_99v99', 't1e_95v95', 't1e_90v90', 't1e_80v80')

t1e_ext_v_ext[1, 1] <- my.power(counts100v100$P.Value)
t1e_ext_v_ext[1, 2] <- my.power(counts99v99$P.Value)
t1e_ext_v_ext[1, 3] <- my.power(counts95v95$P.Value)
t1e_ext_v_ext[1, 4] <- my.power(counts90v90$P.Value)
t1e_ext_v_ext[1, 5] <- my.power(counts80v80$P.Value)

write.csv(t1e_ext_v_ext, paste0(dir, "T1e_NFE_100v100-80v80_maf", maf, ".csv"), quote=F, row.names=F)
