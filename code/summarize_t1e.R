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
# scen = "s1"
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
Ncc = 'cc10k'  #Number of common controls: 'cc5k' or 'cc10k'
int_prune = 100
ext_prune = 100
pruning = "pruneSepRaresim" #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '100v80'
data = 'by_gene'
method = 'skato' #Options: prox, prox_weighted, prox2, iecat, skato, skat, burden
type = 'ext' #Options: int, ext, all

# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/', int_prune, 'v', ext_prune, '/')

### Type 1 error

# t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_adj = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_homo = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_skat_", pruning, "_", data, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T) # for int v ext
# t1e_pruning = read.table(paste0(dir, "T1e_skat_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_gene_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
t1e_pruning = read.table(paste0(dir, "T1e_gene_", method, "_", type, "_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)


# t1e_out = apply(t1e, 2, my.power)
# t1e_adj_out = apply(t1e_adj, 2, my.power)
# t1e_homo_out = apply(t1e_homo, 2, my.power)
# t1e_skat_out = apply(t1e_skat, 2, my.power)
# t1e_og_hap_out = apply(t1e_og_hap, 2, my.power)
t1e_pruning_out = apply(t1e_pruning, 2, my.power)

# write.csv(t(as.data.frame(t1e_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_adj_out)), paste0(dir, "T1e_all_adj_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_homo_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_skat_", pruning, "_", data, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F) # for int v ext
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_skat_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_gene_", method, "_", type, "_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)



dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/') 

prox = read.table(paste0(dir, "T1e_gene_prox_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
prox2 = read.table(paste0(dir, "T1e_gene_prox2_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)

comb_int = cbind(prox_CLEC17A = prox$CLEC17A, prox2_CLEC17A = prox2$CLEC17A, 
                 prox_NDUFB7 = prox$NDUFB7, prox2_NDUFB7 = prox2$NDUFB7)

comb_ext = cbind(prox_ADGRE2 = prox$ADGRE2, prox2_ADGRE2 = prox2$ADGRE2,
                 prox_CLEC17A = prox$CLEC17A, prox2_CLEC17A = prox2$CLEC17A,
                 prox_DNAJB1 = prox$DNAJB1, prox2_DNAJB1 = prox2$DNAJB1,
                 prox_NDUFB7 = prox$NDUFB7, prox2_NDUFB7 = prox2$NDUFB7,
                 prox_PTGER1 = prox$PTGER, prox2_PTGER = prox2$PTGER)


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
