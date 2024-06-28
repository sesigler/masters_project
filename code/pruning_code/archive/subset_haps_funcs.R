############################################################################## 
# This file contains the functions necessary for choosing which individuals
# will be selected for each dataset subset by ancestry as well as functions
# to perform the actual subsetting of the datasets
# Variable legend:
# pop<#>.cols: the vector containing the indices for a specified population
# pop<#>.<DATASET>.size: number of haplotypes required for a particular dataset and population
# scen: which scenario of the pipeline you are running (options: s1, s2)  
# <DATASET>.cols: the column indices that have been selected for a particular dataset
# hap: the haplotype file
# Pop: the ancestry proportion (e.g., AFR, NFE)
##############################################################################


### Subset columns from haplotype file by ancestry for cases
case_cols = function(pop1.cols, pop2.cols, pop1.case.size, pop2.case.size, scen) {
  
  if (scen == 's1') {
    cases.pop1 = sort(sample(x=pop1.cols, size = pop1.case.size, replace = FALSE))
    cases.pop2 = sort(sample(x=pop2.cols, size = pop2.case.size, replace = FALSE))
    
    cases = list(cases.pop1, cases.pop2)
    return(cases)
  }
  else if (scen == 's2') {
    cases.pop1 = sort(sample(x=pop1.cols, size = pop1.case.size, replace = FALSE))
    
    return(cases.pop1)
  }
  else {
    print("Scenario not recognized.")
  }
}

### Subset columns from haplotype file by ancestry for internal controls
int_cols = function(pop1.cols, pop2.cols, pop1.int.size, pop2.int.size, case.cols, scen) {
  
  if (scen == 's1') {
    pop1.case.cols = case.cols[[1]]
    pop2.case.cols = case.cols[[2]]
    
    # Remove elements by value not index
    int.pop1 = sort(sample(x=pop1.cols[! pop1.cols %in% pop1.case.cols], size = pop1.int.size, replace = FALSE))
    int.pop2 = sort(sample(x=pop2.cols[! pop2.cols %in% pop2.case.cols], size = pop2.int.size, replace = FALSE))
    
    ints = list(int.pop1, int.pop2)
    return(ints)
  }
  else if (scen == 's2') {
    int.pop1 = sort(sample(x=pop1.cols[! pop1.cols %in% case.cols], size = pop1.int.size, replace = FALSE))
    
    return(int.pop1)
  }
  else {
    print("Scenario not recognized.")
  }
}

### Subset columns from haplotype file by ancestry for common controls
cc_cols = function(pop1.cols, pop2.cols, pop1.cc.size, pop2.cc.size, case.cols, int.cols, scen) {
  
  if (scen == 's1') {
    pop1.case.cols = case.cols[[1]]
    pop2.case.cols = case.cols[[2]]
    
    pop1.int.cols = int.cols[[1]]
    pop2.int.cols = int.cols[[2]]
    
    cc.pop1 = sort(sample(x=pop1.cols[! pop1.cols %in% c(pop1.case.cols, pop1.int.cols)], size = pop1.cc.size, replace = FALSE))
    cc.pop2 = sort(sample(x=pop2.cols[! pop2.cols %in% c(pop2.case.cols, pop2.int.cols)], size = pop2.cc.size, replace = FALSE))
    
    ccs = list(cc.pop1, cc.pop2)
    return(ccs)
  }
  else if (scen == 's2') {
    cc.pop1 = sort(sample(x=pop1.cols[! pop1.cols %in% c(case.cols, int.cols)], size = pop1.cc.size, replace = FALSE))
    cc.pop2 = sort(sample(x=pop2.cols, size = pop2.cc.size, replace = FALSE))
    
    ccs = list(cc.pop1, cc.pop2)
    return(ccs)
  }
  else {
    print("Scenario not recognized.")
  }
}

### Subset columns from haplotype file by ancestry for references
ref_cols = function(pop1.cols, pop2.cols, pop1.ref.size, pop2.ref.size, case.cols, int.cols, cc.cols, scen) {

  pop1.cc.cols = cc.cols[[1]]
  pop2.cc.cols = cc.cols[[2]]
  
  if (scen == 's1') {
    pop1.case.cols = case.cols[[1]]
    pop2.case.cols = case.cols[[2]]
    
    pop1.int.cols = int.cols[[1]]
    pop2.int.cols = int.cols[[2]]
    
    ref.pop1 = sort(sample(x=pop1.cols[! pop1.cols %in% c(pop1.case.cols, pop1.int.cols, pop1.cc.cols)], size = pop1.ref.size, replace = FALSE))
    ref.pop2 = sort(sample(x=pop2.cols[! pop2.cols %in% c(pop2.case.cols, pop2.int.cols, pop2.cc.cols)], size = pop2.ref.size, replace = FALSE))
    
    refs = list(ref.pop1, ref.pop2)
    return(refs)
  }
  else if (scen == 's2') {
    ref.pop1 = sort(sample(x=pop1.cols[! pop1.cols %in% c(case.cols, int.cols, pop1.cc.cols)], size = pop1.ref.size, replace = FALSE))
    ref.pop2 = sort(sample(x=pop2.cols[! pop2.cols %in% pop2.cc.cols], size = pop2.ref.size, replace = FALSE))
    
    refs = list(ref.pop1, ref.pop2)
    return(refs)
  }
  else {
    print("Scenario not recognized.")
  }
}


# Function to subset the haplotype file based on dataset (cases and internal controls)
sub_hap_scen = function(hap, cols, scen) {
  
  if (scen == 's1') {
    pop1.cols = cols[[1]]
    pop2.cols = cols[[2]]
    
    hap_sub = hap[, c(pop1.cols, pop2.cols)]
    
    return(hap_sub)
  }
  else if (scen == 's2') {
    hap_sub = hap[, cols]
    
    return (hap_sub)
  }
  else {
    print("Scenario not recognized.")
  }
}

# Function to subset the haplotype file based on common controls
sub_hap_cc = function(hap, cols) {
  
    pop1.cols = cols[[1]]
    pop2.cols = cols[[2]]
    
    hap_sub = hap[, c(pop1.cols, pop2.cols)]
    
    return(hap_sub)
}

# Function to subset the haplotype file based on references
sub_hap_refs = function(hap, cols, Pop) {
  
  if (Pop == 'AFR') {
    
    pop.cols = cols[[1]]
    hap_sub = hap[, pop.cols]
    
    return(hap_sub)
  }
  else if (Pop == 'NFE') {
   
    pop.cols = cols[[2]]
    hap_sub = hap[, pop.cols]
    
    return(hap_sub)
  }
  else {
    print("Invalid population specified.")
  }
}
