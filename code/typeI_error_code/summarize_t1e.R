################################################################################
# Summarize type I error or power results 
################################################################################

# define function for calculating power/t1e
my.power <- function(values, alpha=0.05){
  
  # Exclude NA values 
  values2 = values[!is.na(values)]
  
  # See which values are less than the significance level alpha
  sig = which(as.numeric(values2) <= alpha)
  
  # Calculate the number of significant values out of the total number of non NA values
  out = length(sig)/length(values2)
  
  return(out)
}

# Define parameters
calc = 'T1e'
Pops = c('AFR', 'NFE')
admx_props = c(80, 20)
scen = 's2'
sub_scen = 'default'
maf = 0.001 

dir = paste0('/home/math/siglersa/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/', tolower(calc), '/')
# dir = paste0('C:/Users/sagee/Documents/GitHub/masters_project/Data/admixed/', paste(paste(admx_props, Pops, sep = ""), collapse = "_"), '/', scen, '/', sub_scen, '/', tolower(calc), '/')

file_in = paste0(scen, "_", sub_scen, "_maf", maf, ".txt")
file_out = paste0(scen, "_", sub_scen, "_maf", maf, ".csv")


# Get median Neff size
neff_all = read.csv(paste0(dir, "neff", "_", scen, "_", sub_scen, "_maf", maf, ".csv"), header=T)
print(summary(neff_all))

# Read in results for each method
prox_ext = read.table(paste0(dir, calc, "_gene_prox_ext_", file_in), header = T)
prox_ext_adj_Ncc = read.table(paste0(dir, calc, "_gene_prox_ext_adj_Ncc_", file_in), header = T)
prox_ext_adj_Neff = read.table(paste0(dir, calc, "_gene_prox_ext_adj_Neff_", file_in), header = T)

proxW_ext = read.table(paste0(dir, calc, "_gene_prox_weighted_ext_", file_in), header = T)
proxW_ext_adj_Ncc = read.table(paste0(dir, calc, "_gene_prox_weighted_ext_adj_Ncc_", file_in), header = T)
proxW_ext_adj_Neff = read.table(paste0(dir, calc, "_gene_prox_weighted_ext_adj_Neff_", file_in), header = T)

prox2_ext = read.table(paste0(dir, calc, "_gene_prox2_ext_", file_in), header = T)
prox2_ext_adj_Ncc = read.table(paste0(dir, calc, "_gene_prox2_ext_adj_Ncc_", file_in), header = T)
prox2_ext_adj_Neff = read.table(paste0(dir, calc, "_gene_prox2_ext_adj_Neff_", file_in), header = T)

prox2_all = read.table(paste0(dir, calc, "_gene_prox2_all_", file_in), header = T)
prox2_all_adj_Ncc = read.table(paste0(dir, calc, "_gene_prox2_all_adj_Ncc_", file_in), header = T)
prox2_all_adj_Neff = read.table(paste0(dir, calc, "_gene_prox2_all_adj_Neff_", file_in), header = T)

iecat_all = read.table(paste0(dir, calc, "_gene_iecat_all_", file_in), header = T)
iecat_all_adj_Ncc = read.table(paste0(dir, calc, "_gene_iecat_all_adj_Ncc_", file_in), header = T)
iecat_all_adj_Neff = read.table(paste0(dir, calc, "_gene_iecat_all_adj_Neff_", file_in), header = T)

skato_int = read.table(paste0(dir, calc, "_gene_skato_int_", file_in), header = T)
skato_ext = read.table(paste0(dir, calc, "_gene_skato_ext_", file_in), header = T)
skato_all = read.table(paste0(dir, calc, "_gene_skato_all_", file_in), header = T)

skat_int = read.table(paste0(dir, calc, "_gene_skat_int_", file_in), header = T)
skat_ext = read.table(paste0(dir, calc, "_gene_skat_ext_", file_in), header = T)
skat_all = read.table(paste0(dir, calc, "_gene_skat_all_", file_in), header = T)

burden_int = read.table(paste0(dir, calc, "_gene_burden_int_", file_in), header = T)
burden_ext = read.table(paste0(dir, calc, "_gene_burden_ext_", file_in), header = T)
burden_all = read.table(paste0(dir, calc, "_gene_burden_all_", file_in), header = T)

# Calculate t1e or power for each method
prox_ext = apply(prox_ext, 2, my.power)
prox_ext_adj_Ncc = apply(prox_ext_adj_Ncc, 2, my.power)
prox_ext_adj_Neff = apply(prox_ext_adj_Neff, 2, my.power)

proxW_ext = apply(proxW_ext, 2, my.power)
proxW_ext_adj_Ncc = apply(proxW_ext_adj_Ncc, 2, my.power)
proxW_ext_adj_Neff = apply(proxW_ext_adj_Neff, 2, my.power)

prox2_ext = apply(prox2_ext, 2, my.power)
prox2_ext_adj_Ncc = apply(prox2_ext_adj_Ncc, 2, my.power)
prox2_ext_adj_Neff = apply(prox2_ext_adj_Neff, 2, my.power)

prox2_all = apply(prox2_all, 2, my.power)
prox2_all_adj_Ncc = apply(prox2_all_adj_Ncc, 2, my.power)
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

# Combine all results into one dataframe
results = cbind(Gene=names(prox_ext), prox_ext, prox_ext_adj_Ncc, prox_ext_adj_Neff, 
                proxW_ext, proxW_ext_adj_Ncc, proxW_ext_adj_Neff, 
                prox2_ext, prox2_ext_adj_Ncc, prox2_ext_adj_Neff,
                prox2_all, prox2_all_adj_Ncc, prox2_all_adj_Neff,
                iecat_all, iecat_all_adj_Ncc, iecat_all_adj_Neff, 
                skato_int, skato_ext, skato_all,
                skat_int, skat_ext, skat_all,
                burden_int, burden_ext, burden_all)

# Save table of results
write.csv(as.data.frame(results), paste0(dir, calc, "_all_gene_", file_out), quote=F, row.names=F)

