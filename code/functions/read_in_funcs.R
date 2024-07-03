############################################################################## 
# This file contains the functions necessary to read in the legend and haplotype
# files necessary to perform type I error calculations
#
# Variable legend:
#' @param dir: directory where the file is stored
#' @param Pop: Population being simulated, can be a homogeneous population (e.g. AFR)
#'             or an admixed population such as AFR_NFE or LTX
#' @param sim_rep: current simulation replicate in for loop (i.e., i)
#' @param scen: which scenario of the pipeline you are running (options: s1, s2)
#' @param dataset: which dataset the haplotype is coming from 
#                  (options: cases, internal.controls, common.controls, refs)
#' @param p_fun: percent pruned for functional variants
#' @param p_syn: percent pruned for synonymous variants
##############################################################################

# Read in legend file
read_leg = function(dir, Pop, sim_rep) {
  
  leg = read.table(paste0(dir, 'chr19.block37.', Pop, '.sim', sim_rep, '.legend'), header=T, sep='\t')
  leg$row = 1:nrow(leg)
  
  return(leg)
}

# Read in haplotype file (admixed)
read_hap = function(dir, Pop, sim_rep, scen, dataset, p_fun, p_syn) {
  
  hap = fread(paste0(dir, 'chr19.block37.', Pop, '.sim', sim_rep, '.', scen, '.', dataset, '.', p_fun, 'fun.', p_syn, 'syn.haps.gz'))
  hap = as.data.frame(hap)
  
  return(hap)
}

# Read in haplotype file (homogeneous)-does not contain scen parameter
read_hap_homo = function(dir, Pop, sim_rep, dataset, p_fun, p_syn) {
  
  hap = fread(paste0(dir, 'chr19.block37.', Pop, '.sim', sim_rep, '.', dataset, '.', p_fun, 'fun.', p_syn, 'syn.haps.gz'))
  hap = as.data.frame(hap)
  
  return(hap)
}
