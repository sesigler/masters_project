### Jess's typeI_error_gene.R code file

library(dplyr)
library(tidyr)
library(proxecat)
library(SKAT)
library(iECAT)

dir = "/storage/math/projects/compinfo/simulations/output/22K_together_gene/"

# create empty vectors to store the p-values from each replicate
prox.p = prox2.p = iecat.p = skat.p = c()

# loop through the simulation replicates
for (i in 1:100){
  
   # read in files
   leg = read.table(paste0(dir, "chr19.block37.NFE.sim", i, ".legend"), header=T)
   count = read.table(paste0(dir, "chr19.block37.NFE.sim", i, ".count.txt"))
   sample = read.table(paste0(dir, "chr19.block37.NFE.sim", i, ".controls.sample"), header=T) %>% filter(ID_1!=0)
   haps = read.table(gzfile(paste0(dir, "chr19.block37.NFE.sim", i, ".controls.haps.gz")))

   # ensure the phenotype column is numeric
   sample$pheno = as.numeric(as.character(sample$pheno))
      
   # designate cases/controls
   # 0 = internal control, 1 = case, 2 = external control, 3 = external control
   cases = sample(1:nrow(sample), 1000)
   sample$pheno[cases] = 1
   ext.controls = sample(which(sample$pheno==0), 10000)   
   sample$pheno[ext.controls] = 2
   ext.controls2 = sample(which(sample$pheno==0), 10000)   
   sample$pheno[ext.controls2] = 3
   
   # create a genotype matrix
   geno = matrix(0, nrow(haps), ncol(haps)/2)
   
   # sum up the number of alleles in adjacent haplotypes (2 haplotypes per person)
   for (j in 1:(ncol(haps)/2)){
     geno[,j] = haps[,2*j] + haps[,2*j-1]
   }
   geno = as.data.frame(geno)

   # remove the second set of external control samples
   geno2 = geno %>% select(!paste0("V", ext.controls2))

   # calculate minor allele counts (MACs)
   sum = rowSums(geno2)

   # identify/remove the common variants
   N = ncol(geno)-length(ext.controls2) # number of individuals
   maf1 = round(0.01*(2*N)) # minor allele frequency (MAF) of 1%
   common = which(sum>maf1)

   # add annotations and counts to genotype matrix
   geno.prox = cbind(leg[-common, c("gene", "fun")], geno2[-common,])
   
   # convert all 0 counts to NAs (makes next step faster)
   geno.prox.na = na_if(geno.prox, "0")

   # create a dataframe with a line for each variant instead of just counts (necessary for ProxECAT v2)
   data = pivot_longer(geno.prox.na, cols=names(geno2), names_to=c("name"), values_to="count", values_drop_na=T)

   # ensure individuals with two copies of the variant are listed twice
   data2 = data %>% filter(count==2)
   data.all = rbind(data, data2) %>% select(-count)

   # add case/control status and functional group info
   data.all$case = as.factor(ifelse(data.all$name %in% paste0("V", cases), "case", "control"))
   data.all$group = as.factor(ifelse(data.all$name %in% paste0("V", ext.controls), "ext", "int"))
   data.all$fun = as.factor(data.all$fun)

   # count the number of alleles per gene per status (case/control & fun/syn)
   counts = data.all %>% count(gene, case, fun)
   counts.wide = tidyr::pivot_wider(counts, names_from=c(case, fun), values_from=n,
                          values_fill=0, names_sep=".")

   # call proxECAT
   prox = proxecat(counts.wide$case.fun, counts.wide$case.syn, counts.wide$control.fun, counts.wide$control.syn)

   # store proxECAT p-values
   prox.p = rbind(prox.p, prox$p.value)

   # remove external controls
   pheno = sample$pheno[-c(ext.controls, ext.controls2)]

   # create genotype matrix for internal cases/controls & external controls
   geno.skat = cbind(leg[-common, c("gene", "fun")], geno[-common, -c(ext.controls, ext.controls2)]) %>% filter(fun=="fun")
   geno.ext = cbind(leg[-common, c("gene", "fun")], geno[-common, ext.controls]) %>% filter(fun=="fun")

   # add allele counts to external controls matrix
   geno.ext$count = rowSums(geno.ext[,-c(1,2)])

   # null model object
   obj = SKAT_Null_Model(as.numeric(pheno) ~ 1, out_type="D") # D-dichotomous

   # call ProxECATv2/iECAT/SKAT once per gene
   prox2 = iecat = skat = re = c()
   genes = levels(droplevels(as.factor(leg$gene)))
   for(g in 1:length(genes)){ # loop through the genes
  
     # subset the data by gene
     temp = data.all %>% filter(gene==genes[g])

     # return NA if there are no minor alleles in the internal samples or cases
     if (summary(temp$group)[2]==0 | summary(temp$case)[1]==0){

        pvalue = NA
       
     } else {
        
        # fit the ProxECATv2 model
        glm.prox = glm(fun ~ case + group, data=temp, family="binomial") 
  
        # save the p-value for case/control status
        pvalue = summary(glm.prox)$coefficients[2,4]
     }
     
     # genotype matrix for internal cases/controls
     Z = geno.skat %>% filter(gene==genes[g]) %>% select(-gene, -fun)
  
     # MAC matrix for external controls
     tbl = geno.ext %>% filter(gene==genes[g]) %>% select(a0=count) %>% mutate(a1=2*(ncol(geno.ext)-3)-a0)

     # call the iECAT function
     re = iECAT(t(Z), obj, as.matrix(tbl), method="optimal")   

     prox2 = c(prox2, pvalue)
     iecat = c(iecat, re$p.value)
     skat = c(skat, re$p.value.internal)
   }

   # store the ProxECATv2/iECAT/SKAT p-values
   prox2.p = rbind(prox2.p, prox2)
   iecat.p = rbind(iecat.p, iecat)
   skat.p = rbind(skat.p, skat) 

   print(i)
}

dir2 = "/nfs/storage/math/gross-s2/projects/compinfo/simulations/results/"

colnames(prox.p) = colnames(prox2.p) = colnames(iecat.p) = colnames(skat.p) = genes
write.table(prox.p, paste0(dir2, "ProxECAT_TypeIError_", end, ".txt"), quote=F, row.names=F)
write.table(prox2.p, paste0(dir2, "ProxECATv2_TypeIError_", end, ".txt"), quote=F, row.names=F)
write.table(iecat.p, paste0(dir2, "iECAT_TypeIError_", end, ".txt"), quote=F, row.names=F)
write.table(skat.p, paste0(dir2, "SKAT_TypeIError_", end, ".txt"), quote=F, row.names=F)