# Subset from the master legend file of a single population

end = 100 
start = end-99
pop = "EAS"

library(dplyr) 

dir_in = '/storage/math/projects/compinfo/simulations/input/'
dir_out = paste0('/home/math/siglersa/Sim_100k/', pop, '/subset_master/')

# read in master legend file
master = read.table(paste0(dir_in, 'chr19.block37.', pop, '.master.legend'), sep='\t')

colnames(master) = c("position", "id", "a0", "a1", "AC", "prob", "exonic", "gene", "fun")
master$alleles = paste0(master$a0, '/', master$a1)

# subset the file according to the number of alleles at each position
singles = master %>% filter(prob==1)
dups = master %>% filter(prob==0.5)
trips = master %>% filter(prob==".")

for (i in start:end){
  
  # get a list of all the  positions in trips
  trips.po = levels(droplevels(as.factor(trips$position)))
  
  # create a table of transitions/transversions
  trips.po1  =  as.data.frame(matrix(NA, nrow=length(trips.po), ncol=3))
  colnames(trips.po1) = c('position', 'draw', 'TiTv')
  trips.po1$position = trips.po
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
  dups2 = dups %>% group_by(position) %>% sample_n(size=1)
  
  # merge all positions
  master2 = union(singles, union(dups2, trips2)) %>% arrange(position) %>% select(id, position, a0, a1, AC, prob, exonic, gene, fun)
  
  master2$fun = ifelse(master2$fun=="synonymous SNV", "syn", "fun")
  master2$exonic[grepl("exonic", master2$exonic)] = "exonic"
  master2$gene[grepl("ZNF333", master2$gene)] = "ZNF333"
  
  # write output legend file
  write.table(master2, paste0(dir_out, 'chr19.block37.', pop, '.sim', i, '.legend'), 
              row.names=F, col.names=T, quote=F, sep='\t')
}
