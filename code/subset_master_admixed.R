
end = 100 # change back to mysim
start = end-99
set.seed(12345)

# define the admixture proportion
prop_ad = 0.80

library(dplyr)

dir_in = '/storage/math/projects/compinfo/simulations/input/'
dir_out = '/home/math/siglersa/admixed/AFR_NFE_pops/'

# read in the master legend files
master.AFR = read.table(paste0(dir_in, 'chr19.block37.AFR.master.legend'), sep='\t')
master.NFE = read.table(paste0(dir_in, 'chr19.block37.NFE.master.legend'), sep='\t')

# rename the column names
colnames(master.AFR) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
colnames(master.NFE) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
  
# create an alleles column in the format of a0/a1
master.AFR$alleles = paste0(master.AFR$a0, '/', master.AFR$a1)
master.NFE$alleles = paste0(master.NFE$a0, '/', master.NFE$a1)
  
# extract the probabilities at each position for each ancestry
AFR.prob = master.AFR %>% select(position, prob) %>% unique()
NFE.prob = master.NFE %>% select(position, prob) %>% unique()
probs = merge(AFR.prob, NFE.prob, by="position", suffixes=c(".AFR", ".NFE"))
  
# subset the AFR file according to the number of alleles at each position
singles.AFR = master.AFR %>% filter(prob==1)
dups.AFR = master.AFR %>% filter(prob==0.5)
trips.AFR = master.AFR %>% filter(prob==".")
  
# subset the NFE file according to the number of alleles at each position
singles.NFE = master.NFE %>% filter(prob==1)
dups.NFE = master.NFE %>% filter(prob==0.5)
trips.NFE = master.NFE %>% filter(prob==".")
  
# subset the single positions seen in both populations
singles.AFR2 = singles.AFR %>% filter(position %in% singles.NFE$position)
singles.NFE2 = singles.NFE %>% filter(position %in% singles.AFR$position)
  
# keep the single positions with the same alternate allele in both populations
singles.keep = singles.AFR2 %>% filter(alleles==singles.NFE2$alleles)
  
# subset the single positions in just one of the populations
singles.NA.NFE = singles.NFE %>% filter(!(position %in% singles.AFR$position))
singles.NA.AFR = singles.AFR %>% filter(!(position %in% singles.NFE$position))
  
# combine the single positions together
singles.keep2 = rbind(singles.keep, singles.NA.NFE, singles.NA.AFR)
  
# subset the duplicate positions in both populations
dups.po = probs %>% filter((prob.AFR==0.5 & prob.NFE==0.5) | (prob.AFR==0.5 & prob.NFE==".") | (prob.AFR=="." & prob.NFE==0.5)) %>% select(position)
dups.AFR2 = dups.AFR %>% filter(position %in% dups.po$position)
dups.NFE2 = dups.NFE %>% filter(position %in% dups.po$position) %>% filter(!(position %in% dups.AFR2$position))
dups.keep = rbind(dups.AFR2, dups.NFE2)

# subset the triplicate positions in both populations
trips.po = probs %>% filter(prob.AFR=="." & prob.NFE==".") %>% select(position)
trips = trips.AFR %>% filter(position %in% trips.po$position)
    
for (i in start:end){
    
  # create a table of transitions/transversions for the triplicates
  trips.po1  =  as.data.frame(matrix(NA, nrow=nrow(trips.po), ncol=3))
  colnames(trips.po1) = c('position', 'draw', 'TiTv')
  trips.po1$position = trips.po$position
  trips.po1$draw = runif(nrow(trips.po1))
  trips.po1$TiTv[which(trips.po1$draw < 0.7396)] = 'transition'
  trips.po1$TiTv[which(trips.po1$draw >= 0.7396)] = 'transversion'
  #head(trips.po1)
    
  # subset the transitions
  ti = trips.po1 %>% filter(TiTv  == 'transition')
  ann.ti = trips %>% filter(position %in% ti$position, alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C')
    
  # remove all of the transitions
  ann.tv = trips %>% filter(!(position %in% ti$position), !(alleles=='A/G' | alleles=='G/A' | alleles=='C/T' | alleles=='T/C'))
    
  # randomly pick an allele from the transversions
  ann.tv2 = ann.tv %>% group_by(position) %>% sample_n(size=1)
    
  # merge transitions and transversions
  trips2 = union(ann.ti, ann.tv2)
    
  # randomly pick an allele from the duplicates
  dups2 = dups.keep %>% group_by(position) %>% sample_n(size=1)
    
  # subset the single positions with different alternate alleles in the populations
  singles.dif = singles.AFR2 %>% filter(alleles!=singles.NFE2$alleles)
    
  # pick an allele for the single positions with different alternate alleles based on the admixed proportion
  singles.dif.AFR = singles.dif %>% sample_n(size=round(prop_ad*nrow(singles.dif)))
  singles.dif.NFE = singles.NFE2 %>% filter(alleles!=singles.AFR2$alleles) %>% filter(!(position %in% singles.dif.AFR$position))
    
  # combine the single positions together
  singles2 = rbind(singles.keep2, singles.dif.AFR, singles.dif.NFE)
    
  # merge all positions
  master2 = union(singles2, union(dups2, trips2)) %>% arrange(position) %>% 
    select(id, position, a0, a1, AC, prob, exonic, gene, fun)

  master2$fun = ifelse(master2$fun=="synonymous SNV", "syn", "fun")
  master2$exonic[grepl("exonic", master2$exonic)] = "exonic"
  master2$gene[grepl("ZNF333", master2$gene)] = "ZNF333"
    
  # write output legend file
  out.name = paste0(dir_out, "chr19.block37.AFR_NFE.sim", i, ".copy.legend")
  write.table(master2, out.name, row.names=F, col.names=T, quote=F, sep='\t')

  print(i)
  #print(nrow(master2))
}