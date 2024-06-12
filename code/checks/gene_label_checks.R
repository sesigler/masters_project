library(dplyr)
library(tidyr)
library(data.table)
library(proxecat)
library(SKAT)
library(iECAT)

source("/home/math/siglersa/mastersProject/Input/read_in_funcs.R")
source("/home/math/siglersa/mastersProject/Input/general_data_manip.R")
source("/home/math/siglersa/mastersProject/Input/methods_funcs.R")

# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/read_in_funcs.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/general_data_manip.R")
# source("C:/Users/sagee/Documents/GitHub/masters_project/code/typeI_error_code/methods_funcs.R")

Pop = 'NFE'
pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '100v80'
data = 'by_gene'
p_case = 100
p_case_fun = p_case_syn = p_int_fun = p_int_syn = int_prune = 100
p_cc_fun = p_cc_syn = ext_prune = 80
Ncase = Nint = 5000
Ncc = 10000 #Number of common controls: 5000 or 10000 
maf = 0.001 #MAF: 0.001 (0.1%) or 0.01 (1%)
# scen = 's1'

# dir = '/storage/math/projects/compinfo/simulations/'
# setwd(paste0(dir, 'output/20K_', Pop, '/'))

dir_leg = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/attempt2_combine_MACbins_legFiles_differ/')
# dir_in = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/')
dir_out = paste0('/home/math/siglersa/mastersProject/Output/')

# dir_leg = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_in = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
# dir_out = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/output/')

# create empty vectors to store the p-values from each replicate
# prox.p = prox.int.p = prox.all.p = c()
# prox2.p = prox2.int.p = prox2.all.p = c()
# iecat.p = skat.int.p = skat.ext.p = skat.all.p = c()
# prox_genes_p = prox2_genes_p = iecat_genes_p = skat_int_genes_p = skat_ext_genes_p = c()

gene_labels = c() #Burden


# loop through the simulation replicates
set.seed(1) 
for (i in 1:100){
  
  leg = read.table(paste0(dir_leg, 'chr19.block37.', Pop, '.sim', i, '.', p_case, 'fun.', p_case, 'syn.legend'), header=T, sep='\t') #RAREsim v2.1.1 pruning only
  leg$row = 1:nrow(leg)
  
  genes = levels(droplevels(as.factor(leg$gene)))
  gene_labels = c(gene_labels, list(genes))
  
  print(i)
  
}

max_length <- max(sapply(gene_labels, length))

for (i in 1:length(gene_labels)) {
  current_length <- length(gene_labels[[i]])
  if (current_length < max_length) {
    gene_labels[[i]] <- c(gene_labels[[i]], rep(NA, max_length - current_length))
  }
}

results <- do.call(rbind.data.frame, gene_labels)
write.table(results, paste0(dir_out, "gene_labels.txt"), quote=F, row.names=F, col.names=F)
