# Get the downsampling info for the AMR population

# load libraries
library(dplyr)
library(tidyr)
library(data.table)

dir_in = '/storage/math/projects/compinfo/simulations/input/'
dir_out = '/home/math/siglersa/IAM_haps/gnomad/'
# dir_in = dir_out = 'C:/Users/sagee/Documents/HendricksLab/IAM_hap_data/'

# Read in byGene downsampling file
byGene = read.table(paste0(dir_in, "Block37_fun_syn_num_var_byGene.txt"), header = T)

# Sum obs_fun and obs_syn over all genes in each downsample level for each pop
# Add block column
byGene2 = byGene %>% group_by(pop, downsampling) %>% 
  summarize(obs_syn2 = sum(obs_syn), obs_fun2 = sum(obs_fun)) %>%
  mutate(block = "37")

# Output data in same order as Block37_fun_syn_num_var.txt
out = byGene2 %>% select(block, pop, downsample=downsampling, obs_fun=obs_fun2, obs_syn=obs_syn2)

# Save downsampling file with AMR data
write.table(out, paste0(dir_out, "Block37_fun_syn_num_var.txt"), quote=F, row.names=F, col.names=T)

