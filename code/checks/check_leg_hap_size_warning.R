### File to check the sizes of:
# Megan's chr19.block37.NFE.sim${rep}.controls.haps.gz
# Megan's chr19.block37.NFE.sim${rep}.controls.haps.sm
# Megan's chr19.block37.NFE.sim${rep}.copy.legend
# My chr19.block37.NFE.sim${rep}.all.${pcase}fun.${pcase}syn.haps.gz
# My chr19.block37.NFE.sim${rep}.all.${pcase}fun.${pcase}syn.haps.sm
# My chr19.block37.NFE.sim${rep}.${pcase}fun.${pcase}syn.legend

#Added comment based on Adelle's peer review 
library(data.table)
library(dplyr)

Pop = 'NFE'
p_case = 100
# p_exp = 100
# p_conf = 80
# Nsim = 20000 
pruning = 'pruneSepRaresim' #Options: pruneSeparately, pruneSequentially, pruneTogether, pruneSepRaresim, pruneSepR
folder = '100v80'
# int_prune = 100
# ext_prune = 80

# Added another comment based on Adelle's peer review suggestion
meg_dir = '/storage/math/projects/RAREsim/Cases/Sim_20k/NFE/data/'
my_dir = paste0('/home/math/siglersa/mastersProject/20K_NFE/', pruning, '/', folder, '/attempt2_combine_MACbins_legFiles_differ/')

meg_dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')
my_dir = paste0('C:/Users/sagee/Documents/HendricksLab/mastersProject/input/', pruning, '/', folder, '/')

comp <- data.frame(matrix(ncol = 7, nrow = 100))
colnames(comp) <- c('mn_hap_gz_length', 'mn_hap_sm_length', 'mn_leg_length',
                    'my_hap_gz_length', 'my_hap_sm_length', 'my_leg_length',
                    'rep')

j <- 7
for (j in 1:100) {
  
  # read in Megan's files
  meg_hap_gz = fread(paste0(meg_dir, 'chr19.block37.', Pop, '.sim', j, '.controls.haps.gz'))
  # meg_hap_sm = fread(paste0(meg_dir, 'chr19.block37.', Pop, '.sim', j, '.controls.haps.sm'))
  meg_leg = read.table(paste0(meg_dir, 'chr19.block37.', Pop, '.sim', j, '.copy.legend'), header=T, sep='\t') 
  
  # read in my files
  my_hap_gz = fread(paste0(my_dir, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.', p_case, 'syn.haps.gz'))
  # my_hap_sm = fread(paste0(my_dir, 'chr19.block37.', Pop, '.sim', j, '.all.', p_case, 'fun.', p_case, 'syn.haps.sm'))
  my_leg = read.table(paste0(my_dir, 'chr19.block37.', Pop, '.sim', j, '.', p_case, 'fun.', p_case, 'syn.legend'), header=T, sep='\t')
  
  comp[j, 1] <- nrow(meg_hap_gz)
  comp[j, 2] <- nrow(meg_hap_sm)
  comp[j, 3] <- nrow(meg_leg)
  
  comp[j, 4] <- nrow(my_hap_gz)
  comp[j, 5] <- nrow(my_hap_sm)
  comp[j, 6] <- nrow(my_leg)
  
  comp[j, 7] <- j
}
