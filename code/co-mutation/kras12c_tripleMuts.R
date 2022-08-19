library(here)
library(genie.6.1)
library(tidyverse)
library(dplyr)

data(mutation_data)
data(ras_mut_allele)

out.dir = here("output", "co-mutation")
#####

#functions
gene_samples = function(gene, muts){
  return(dim(subset(muts, hugo_symbol == gene))[1])
}



gene_intersect = function(gene1, gene2, muts){
  
  print(paste(gene1, gene2,sep =  "-"))
  gene1.samples = filter(muts, hugo_symbol == gene1)$sample_id
  gene2.samples = filter(muts, hugo_symbol == gene2)$sample_id
  intersect.samples = intersect(gene1.samples, gene2.samples)
  return(length(intersect.samples)) }



#######
kras_12c = ras_mut_allele["KRAS_12C", ]
kras_12c_ids = names(which(kras_12c == 1))

kras_12c_muts = subset(mutation_data, is_matched) %>%
  subset(., sample_id %in% kras_12c_ids)

kras_12c_cancer = as.data.frame(table(kras_12c_muts$cancer)) %>%
  arrange(desc(Freq))


master = c()
for(cancer_type in unique(kras_12c_muts$cancer)){
  print(paste0("Processing " , cancer_type))
  
  muts = subset(kras_12c_muts, cancer == cancer_type)
  gene.muts = data.frame(gene = unique(muts$hugo_symbol))  %>%
    slice(-which(gene == "KRAS")) %>%
    rowwise()  %>%
    mutate(n = gene_samples(gene, muts)) %>%
    filter(n>=5)
  
  if(dim(gene.muts)[1]>1){       
    gene.comb = do.call(expand.grid, rep(list(unique(gene.muts$gene)), 2)) %>%
      dplyr::rename(gene1 = Var1, gene2 = Var2) %>%
      slice(- which(gene2 == gene1))  
    gene.comb <- data.frame(t(apply(gene.comb,1,sort)))
    gene.comb <- gene.comb[!duplicated(gene.comb),]
    
    gene.comb <- gene.comb %>% 
      dplyr::rename(gene1 = X1, gene2 = X2) %>%        
      rowwise()  %>%
      mutate(intersect.n = gene_intersect(gene1,gene2,muts)) 
    # filter(intersect.n >=5)
    master[[cancer_type]] = gene.comb
  }
  
} 

kras12c_tripleMuts = c()
for(ind in c(1:length(master))){
  total.n = dim(subset(kras_12c_muts, cancer == names(master)[[ind]]))[1]
  tmp.df = master[[ind]]
  newDf <- data.frame(t(apply(tmp.df,1,sort)))
  newDf <- newDf[!duplicated(newDf),]
  newDf = newDf %>%
    arrange(desc(intersect.n)) %>%
    mutate(gene = paste0(gene1, "_", gene2)) %>%
    mutate(intersect.percent = round( (100*as.numeric(intersect.n))/total.n, 2) ) %>%
    select(-c(gene1, gene2))
  kras12c_tripleMuts[[names(master)[ind]]] = newDf
  
}

save(kras12c_tripleMuts, file = file.path(out.dir, "kras12c_tripleMuts.rda"))

