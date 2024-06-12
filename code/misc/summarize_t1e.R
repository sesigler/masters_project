#####
# Summarize type I error and/or power results for a particular iteration of the code
#####

# define function for calculating power
my.power <- function(values, alpha=0.05){
  values2 = values[!is.na(values)]
  sig = which(as.numeric(values2) <= alpha)
  out = length(sig)/length(values2)
  return(out)
}

library(stringr) # for str_to_title
# Pop = 'AFR'
# Nsim = '20K'
# pruning = "pruneSepRaresim" #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
# admx = '80-20' #AFR/NFE
calc = 'power'
Pop1 = 'AFR'
Pop2 = 'NFE'
admx_pop1 = 80
admx_pop2 = 20
Nsim = '42k'
scen = 's2'
folder = '160v100v80'
int_prune = 100
ext_prune = 80
Ncase = Nic = 5000
Ncc = 10000 
Nref_pop1 = 10000
Nref_pop2 = 10000
maf = 0.001 
sim_params = paste0('Ncase', Ncase, '_Nic', Nic, '_Ncc', Ncc, '_', Pop1, 'ref', Nref_pop1, '_', Pop2, 'ref', Nref_pop2)
# data = 'by_gene'
# method = 'skato' #Options: prox, prox_weighted, prox2, iecat, skato, skat, burden
# type = 'ext' #Options: int, ext, all


# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', Nsim, '_', Pop, '/', data, '/', folder, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', Pop1, '_', Pop2, '_pops/Sim_', Nsim, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params, '/prox_gene_adj_', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/')
dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', admx_pop1, Pop1, '_', admx_pop2, Pop2, '/Sim_', Nsim, '/', sim_params, '/', scen, '_', folder, '_', int_prune, 'v', ext_prune, '/', calc, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/', pruning, '/', data, '/', folder, '/', int_prune, 'v', ext_prune, '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/checks/')

file_in = paste0(int_prune, "_v_", ext_prune, "_", Pop1, "_", Pop2, "_", scen, "_maf", maf, ".txt")
file_out = paste0(admx_pop1, Pop1, "_", admx_pop2, Pop2, "_", int_prune, "v", ext_prune, "_maf", maf, "_", scen, ".csv")

### Type 1 error

# t1e = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_adj = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", scen, "_adj_", Pop1, '-', Pop2, "_maf", maf, ".txt"), header = T)
# t1e_homo = read.table(paste0(dir, "T1e_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_skat_", pruning, "_", data, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T) # for int v ext
# t1e_pruning = read.table(paste0(dir, "T1e_skat_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_gene_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# t1e_pruning = read.table(paste0(dir, "T1e_gene_", method, "_", type, "_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)

# Get median Neff size
neff_all = read.csv(paste0(dir, "neff.csv"), header=T)
summary(neff_all)

# By gene
prox_int = read.table(paste0(dir, str_to_title(calc), "_gene_prox_int_", file_in), header = T)
prox_ext = read.table(paste0(dir, str_to_title(calc), "_gene_prox_ext_", file_in), header = T)
prox_ext_adj_Ncc = read.table(paste0(dir, str_to_title(calc), "_gene_prox_ext_adj_Ncc_", file_in), header = T)
prox_ext_adj_Neff = read.table(paste0(dir, str_to_title(calc), "_gene_prox_ext_adj_Neff_", file_in), header = T)
# prox_ext_var_adj_Ncc = read.table(paste0(dir, "T1e_gene_prox_ext_var_adj_Ncc_", file_in), header = T)
# prox_ext_var_adj_Neff = read.table(paste0(dir, "T1e_gene_prox_ext_var_adj_Neff_", file_in), header = T)
# prox_ext_gene_adj_Ncc = read.table(paste0(dir, "T1e_gene_prox_ext_gene_adj_Ncc_", file_in), header = T)
# prox_ext_gene_adj_Neff = read.table(paste0(dir, "T1e_gene_prox_ext_gene_adj_Neff_", file_in), header = T)
# prox_ext_adj = read.table(paste0(dir, "T1e_gene_prox_ext_adj_", file_in), header = T)

proxW_int = read.table(paste0(dir, str_to_title(calc), "_gene_prox_weighted_int_", file_in), header = T)
proxW_ext = read.table(paste0(dir, str_to_title(calc), "_gene_prox_weighted_ext_", file_in), header = T)
proxW_ext_adj_Ncc = read.table(paste0(dir, str_to_title(calc), "_gene_prox_weighted_ext_adj_Ncc_", file_in), header = T)
proxW_ext_adj_Neff = read.table(paste0(dir, str_to_title(calc), "_gene_prox_weighted_ext_adj_Neff_", file_in), header = T)
# proxW_ext_var_adj_Ncc = read.table(paste0(dir, "T1e_gene_prox_weighted_ext_var_adj_Ncc_", file_in), header = T)
# proxW_ext_var_adj_Neff = read.table(paste0(dir, "T1e_gene_prox_weighted_ext_var_adj_Neff_", file_in), header = T)
# proxW_ext_gene_adj_Ncc = read.table(paste0(dir, "T1e_gene_prox_weighted_ext_gene_adj_Ncc_", file_in), header = T)
# proxW_ext_gene_adj_Neff = read.table(paste0(dir, "T1e_gene_prox_weighted_ext_gene_adj_Neff_", file_in), header = T)
# proxW_ext_adj = read.table(paste0(dir, "T1e_gene_prox_weighted_ext_adj_", file_in), header = T)

prox2_int = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_int_", file_in), header = T)
prox2_ext = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_ext_", file_in), header = T)
prox2_all = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_all_", file_in), header = T)
prox2_ext_adj_Ncc = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_ext_adj_Ncc_", file_in), header = T)
prox2_ext_adj_Neff = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_ext_adj_Neff_", file_in), header = T)
prox2_all_adj_Ncc = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_all_adj_Ncc_", file_in), header = T)
prox2_all_adj_Neff = read.table(paste0(dir, str_to_title(calc), "_gene_prox2_all_adj_Neff_", file_in), header = T)

iecat_all = read.table(paste0(dir, str_to_title(calc), "_gene_iecat_all_", file_in), header = T)
iecat_all_adj_Ncc = read.table(paste0(dir, str_to_title(calc), "_gene_iecat_all_adj_Ncc_", file_in), header = T)
iecat_all_adj_Neff = read.table(paste0(dir, str_to_title(calc), "_gene_iecat_all_adj_Neff_", file_in), header = T)

skato_int = read.table(paste0(dir, str_to_title(calc), "_gene_skato_int_", file_in), header = T)
skato_ext = read.table(paste0(dir, str_to_title(calc), "_gene_skato_ext_", file_in), header = T)
skato_all = read.table(paste0(dir, str_to_title(calc), "_gene_skato_all_", file_in), header = T)

skat_int = read.table(paste0(dir, str_to_title(calc), "_gene_skat_int_", file_in), header = T)
skat_ext = read.table(paste0(dir, str_to_title(calc), "_gene_skat_ext_", file_in), header = T)
skat_all = read.table(paste0(dir, str_to_title(calc), "_gene_skat_all_", file_in), header = T)

burden_int = read.table(paste0(dir, str_to_title(calc), "_gene_burden_int_", file_in), header = T)
burden_ext = read.table(paste0(dir, str_to_title(calc), "_gene_burden_ext_", file_in), header = T)
burden_all = read.table(paste0(dir, str_to_title(calc), "_gene_burden_all_", file_in), header = T)


# t1e_out = apply(t1e, 2, my.power)
# t1e_adj_out = apply(t1e_adj, 2, my.power)
# t1e_homo_out = apply(t1e_homo, 2, my.power)
# t1e_skat_out = apply(t1e_skat, 2, my.power)
# t1e_og_hap_out = apply(t1e_og_hap, 2, my.power)
# t1e_pruning_out = apply(t1e_pruning, 2, my.power)

# By gene
prox_int = apply(prox_int, 2, my.power)
prox_ext = apply(prox_ext, 2, my.power)
prox_ext_adj_Ncc = apply(prox_ext_adj_Ncc, 2, my.power)
prox_ext_adj_Neff = apply(prox_ext_adj_Neff, 2, my.power)
# prox_ext_var_adj_Ncc = apply(prox_ext_var_adj_Ncc, 2, my.power)
# prox_ext_var_adj_Neff = apply(prox_ext_var_adj_Neff, 2, my.power)
# prox_ext_gene_adj_Ncc = apply(prox_ext_gene_adj_Ncc, 2, my.power)
# prox_ext_gene_adj_Neff = apply(prox_ext_gene_adj_Neff, 2, my.power)
# prox_ext_adj = apply(prox_ext_adj, 2, my.power)

proxW_int = apply(proxW_int, 2, my.power)
proxW_ext = apply(proxW_ext, 2, my.power)
proxW_ext_adj_Ncc = apply(proxW_ext_adj_Ncc, 2, my.power)
proxW_ext_adj_Neff = apply(proxW_ext_adj_Neff, 2, my.power)
# proxW_ext_var_adj_Ncc = apply(proxW_ext_var_adj_Ncc, 2, my.power)
# proxW_ext_var_adj_Neff = apply(proxW_ext_var_adj_Neff, 2, my.power)
# proxW_ext_gene_adj_Ncc = apply(proxW_ext_gene_adj_Ncc, 2, my.power)
# proxW_ext_gene_adj_Neff = apply(proxW_ext_gene_adj_Neff, 2, my.power)
# proxW_ext_adj = apply(proxW_ext_adj, 2, my.power)

prox2_int = apply(prox2_int, 2, my.power)
prox2_ext = apply(prox2_ext, 2, my.power)
prox2_all = apply(prox2_all, 2, my.power)
prox2_ext_adj_Ncc = apply(prox2_ext_adj_Ncc, 2, my.power)
prox2_all_adj_Ncc = apply(prox2_all_adj_Ncc, 2, my.power)
prox2_ext_adj_Neff = apply(prox2_ext_adj_Neff, 2, my.power)
prox2_all_adj_Neff = apply(prox2_all_adj_Neff, 2, my.power)

iecat_all = apply(iecat_all, 2, my.power)
iecat_all_adj_Ncc = apply(iecat_all_adj_Ncc, 2, my.power)
iecat_all_adj_Neff = apply(iecat_all_adj_Neff, 2, my.power)

skato_int = apply(skato_int, 2, my.power)
skato_ext = apply(skato_ext, 2, my.power)
skato_all = apply(skato_all, 2, my.power)

skat_int = apply(skat_int, 2, my.power)
skat_ext = apply(skat_ext, 2, my.power)
skat_all = apply(skat_all, 2, my.power)

burden_int = apply(burden_int, 2, my.power)
burden_ext = apply(burden_ext, 2, my.power)
burden_all = apply(burden_all, 2, my.power)


results = cbind(Gene=names(prox_int), prox_int, prox_ext, prox_ext_adj_Ncc, prox_ext_adj_Neff, 
                proxW_int, proxW_ext, proxW_ext_adj_Ncc, proxW_ext_adj_Neff, 
                prox2_int, prox2_ext, prox2_all, prox2_ext_adj_Ncc, prox2_ext_adj_Neff, prox2_all_adj_Ncc, prox2_all_adj_Neff,
                iecat_all, iecat_all_adj_Ncc, iecat_all_adj_Neff, 
                skato_int, skato_ext, skato_all,
                skat_int, skat_ext, skat_all,
                burden_int, burden_ext, burden_all)

results = cbind(Gene=names(prox_int), prox_int, prox_ext, prox_ext_adj_Ncc, prox_ext_adj_Neff, 
                proxW_int, proxW_ext, proxW_ext_adj_Ncc, proxW_ext_adj_Neff, 
                prox2_int, prox2_ext, prox2_all, prox2_ext_adj_Ncc, prox2_ext_adj_Neff, prox2_all_adj_Ncc, prox2_all_adj_Neff,
                iecat_all, iecat_all_adj_Ncc, iecat_all_adj_Neff)

results = cbind(Gene=names(prox_int), prox_int, prox_ext, prox_ext_adj,
                proxW_int, proxW_ext, proxW_ext_adj,
                prox2_int, prox2_ext, prox2_all, prox2_ext_adj, prox2_all_adj,
                iecat_all, iecat_all_adj, 
                skato_int, skato_ext, skato_all,
                skat_int, skat_ext, skat_all,
                burden_int, burden_ext, burden_all)

results = cbind(Gene=names(prox_int), 
                prox_int, prox_ext, prox_ext_var_adj_Ncc, prox_ext_var_adj_Neff, prox_ext_gene_adj_Ncc, prox_ext_gene_adj_Neff,
                proxW_int, proxW_ext, proxW_ext_var_adj_Ncc, proxW_ext_var_adj_Neff, proxW_ext_gene_adj_Ncc, proxW_ext_gene_adj_Neff)

# results = cbind(Gene=names(prox_ext), prox_ext, prox_ext_adj, proxW_ext, proxW_ext_adj)
# 
# results = cbind(Gene=names(prox_int), prox_int, prox_ext, proxW_int, proxW_ext,
#                 prox2_int, prox2_ext, prox2_all,
#                 iecat_all, skato_int, skato_ext, skato_all,
#                 skat_int, skat_ext, skat_all,
#                 burden_int, burden_ext, burden_all)

# write.csv(t(as.data.frame(t1e_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_adj_out)), paste0(dir, "T1e_all_adj_", int_prune, "_v_", ext_prune, "_", scen, "_", Pop1, '-', Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_homo_out)), paste0(dir, "T1e_all_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_skat_", pruning, "_", data, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F) # for int v ext
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_skat_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(t(as.data.frame(t1e_pruning_out)), paste0(dir, "T1e_all_gene_", method, "_", type, "_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
# write.csv(as.data.frame(results), paste0(dir, "T1e_power_all_gene_", int_prune, "_v_", ext_prune, "_", Pop, "_", Ncc, "_maf", maf, ".csv"), quote=F, row.names=F)
write.csv(as.data.frame(results), paste0(dir, str_to_title(calc), "_all_gene_", file_out), quote=F, row.names=F)
# write.csv(as.data.frame(results), paste0(dir, "T1e_all_prox_gene_adj_", file_out), quote=F, row.names=F)
# write.csv(as.data.frame(results), paste0(dir, "T1e_all_gene_prox_Sim_42k_", file_out), quote=F, row.names=F)



# dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/') 
# 
# prox = read.table(paste0(dir, "T1e_gene_prox_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# prox2 = read.table(paste0(dir, "T1e_gene_prox2_int_", pruning, "_", int_prune, "_v_", ext_prune, "_", Pop2, "_maf", maf, ".txt"), header = T)
# 
# comb_int = cbind(prox_CLEC17A = prox$CLEC17A, prox2_CLEC17A = prox2$CLEC17A, 
#                  prox_NDUFB7 = prox$NDUFB7, prox2_NDUFB7 = prox2$NDUFB7)
# 
# comb_ext = cbind(prox_ADGRE2 = prox$ADGRE2, prox2_ADGRE2 = prox2$ADGRE2,
#                  prox_CLEC17A = prox$CLEC17A, prox2_CLEC17A = prox2$CLEC17A,
#                  prox_DNAJB1 = prox$DNAJB1, prox2_DNAJB1 = prox2$DNAJB1,
#                  prox_NDUFB7 = prox$NDUFB7, prox2_NDUFB7 = prox2$NDUFB7,
#                  prox_PTGER1 = prox$PTGER, prox2_PTGER = prox2$PTGER)
# 
# 
# # Recalculate p-values proxECAT only data
# dir = 'C:/Users/sagee/Documents/GitHub/masters_project/Data/'
# counts99 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_99.csv'), header = T, sep = ',')
# counts95 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_95.csv'), header = T, sep = ',')
# counts90 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_90.csv'), header = T, sep = ',')
# counts80 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_80.csv'), header = T, sep = ',')
# pow <- data.frame(matrix(nrow = 1, ncol = 4))
# colnames(pow) <- c('t1e_99', 't1e_95', 't1e_90', 't1e_80')
# pow[1, 1] <- my.power(counts99$P.Value)
# pow[1, 2] <- my.power(counts95$P.Value)
# pow[1, 3] <- my.power(counts90$P.Value)
# pow[1, 4] <- my.power(counts80$P.Value)
# 
# write.csv(pow, paste0(dir, "T1e_NFE_99-80_maf", maf, ".csv"), quote=F, row.names=F)
# 
# counts100v100 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_100_v_100.csv'), header = T, sep = ',')
# counts99v99 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_99_v_99.csv'), header = T, sep = ',')
# counts95v95 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_95_v_95.csv'), header = T, sep = ',')
# counts90v90 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_90_v_90.csv'), header = T, sep = ',')
# counts80v80 = read.table(paste0(dir, 'proxECAT_counts_expanded_', Pop2, '_80_v_80.csv'), header = T, sep = ',')
# 
# t1e_ext_v_ext <- data.frame(matrix(nrow = 1, ncol = 5))
# colnames(t1e_ext_v_ext) <- c('t1e_100v100', 't1e_99v99', 't1e_95v95', 't1e_90v90', 't1e_80v80')
# 
# t1e_ext_v_ext[1, 1] <- my.power(counts100v100$P.Value)
# t1e_ext_v_ext[1, 2] <- my.power(counts99v99$P.Value)
# t1e_ext_v_ext[1, 3] <- my.power(counts95v95$P.Value)
# t1e_ext_v_ext[1, 4] <- my.power(counts90v90$P.Value)
# t1e_ext_v_ext[1, 5] <- my.power(counts80v80$P.Value)
# 
# write.csv(t1e_ext_v_ext, paste0(dir, "T1e_NFE_100v100-80v80_maf", maf, ".csv"), quote=F, row.names=F)
