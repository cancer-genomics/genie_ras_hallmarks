library(devtools)
library(magrittr)
library(tidyverse)
library(broom)
library(contingency.table)
library(scales)
library(RColorBrewer)
library(viridis)
library(here)
source(here("code", "co-mutation", "stan_funs.R"))
#data(known_drivers, package="rasfuns")

extdir <- system.file("extdata", package="genie.6.1")
load(file.path(extdir, "known_drivers.rda"))

outdir <- here("output", "co-mutation", "summarize_allele_models.R")
fs::dir_create(outdir)
## co-mutation tables
nonras_genes <- function(ct){
  ct$data[[1]] %>%
    pull(nonras) %>%
    unique()
}

metadata <- tibble(script=c("mutations_allele.R",
                            "inactivating_allele.R",
                            "driver_allele.R",
                            "mutations_allele_mutsig.R",
                            "inactivating_allele_mutsig.R",
                             "driver_allele_mutsig.R",
			     "mutations_allele_tmb.R",
                            "inactivating_allele_tmb.R", 
			     "driver_allele_tmb.R", 
			    "mutations_allele_gp.R",
                            "driver_allele_gp.R",
                            "mutations_allele_gp_mutsig.R",
                             "driver_allele_gp_mutsig.R",
                             "mutations_allele_gp_tmb.R",
                             "driver_allele_gp_tmb.R", 
			      "driver_cosmic_allele.R", 
                              "driver_cosmic_allele_gp.R",
                             "driver_cosmic_allele_mutsig.R",
                             "driver_cosmic_allele_gp_mutsig.R",
                             "driver_cosmic_allele_tmb.R",
                             "driver_cosmic_allele_gp_tmb.R"),

                   category=c("Any mutation",
                              "Inactivating mutation",
			      "Driver mutation"	,
                              "Signature",
                              "Inactivating Signature",
			      "Driver Signature",
                              "TMB",
                              "Inactivating TMB", 
			      "Driver TMB", 
			      "Any mutation (gene_pathways)", 
                              "Driver mutation (gene_pathways)",
                              "Signature (gene_pathways)", 
			       "Driver Signature(gene_pathways)", 
		        	 "TMB (gene_pathways)", 
                             "Driver TMB (gene_pathways)", 
                             "Cosmic Driver mutation" , 
                             "Cosmic Driver mutation (gene_pathways)", 
                              "Cosmic Driver Signature",
                              "Cosmic Driver Signature (gene_pathways)",
                              "Cosmic Driver TMB",
                              "Cosmic Driver TMB (gene_pathways)"))

metadata$object <- list(mutation_allele_ct,
                        inactivating_allele_ct,
                        driver_allele_ct, 
			mutation_allele_msig,
                        inactivating_allele_msig,
     			driver_allele_msig,  
                        mutation_allele_tmb,
                        inactivating_allele_tmb, 
			driver_allele_tmb, 
			mutation_allele_gp_ct, 
			driver_allele_gp_ct,
			mutation_allele_gp_msig, 
			driver_allele_gp_msig,  
			mutation_allele_gp_tmb,
			driver_allele_gp_tmb, 
                        driver_cosmic_allele_ct,
                        driver_cosmic_allele_gp_ct, 
                        driver_cosmic_allele_msig,
                        driver_cosmic_allele_gp_msig,
                        driver_cosmic_allele_tmb,
                        driver_cosmic_allele_gp_tmb )

rename <- dplyr::rename

mt <- metadata[c(1,2,3,10,11, 16, 17), ] %>%
  pmap_dfr(harmonize_stats, drivers=known_drivers) %>%
  unnest(stan)

mt2 <- mt %>%
  group_by(group, ras, nonras) %>%
  nest()

mt3 <- metadata[-c(1,2,3,10,11, 16, 17), ] %>%
  pmap_dfr(harmonize_stats2, drivers=known_drivers) %>%
  unnest(stan) %>%
  mutate(category=case_when(category=="TMB" ~ "Any mutation",				
                            category=="Signature"~"Any mutation",
			    category=="TMB (gene_pathways)" ~ "Any mutation (gene_pathways)",		
			    category=="Signature (gene_pathways)"~"Any mutation (gene_pathways)",	
                            category== "Inactivating TMB"~"Inactivating mutation",
                            category== "Inactivating Signature"~"Inactivating mutation",
			category== "Driver TMB (gene_pathways)"~"Driver mutation (gene_pathways)",
                            category== "Driver Signature(gene_pathways)"~"Driver mutation (gene_pathways)",
                             category== "Driver TMB"~"Driver mutation",
                            category== "Driver Signature"~"Driver mutation",
			  category== "Cosmic Driver TMB"~"Cosmic Driver mutation",
                            category== "Cosmic Driver Signature"~"Cosmic Driver mutation",
                        category== "Cosmic Driver TMB (gene_pathways)"~"Cosmic Driver mutation (gene_pathways)",
                            category== "Cosmic Driver.signature(gene_pathways)"~"Cosmic Driver mutation (gene_pathways)",
					TRUE~category))

mt4 <- mt3 %>%
  group_by(group, ras, nonras) %>%
  nest()
combined <- bind_rows(mt2, mt4)
saveRDS(combined, file=here(outdir, "results.rds"))

multicancer <- combined %>%
  ungroup() %>%
  filter(group=="overall") %>%
  rename(cancer_group=group) %>%
  mutate(cancer_group="multi") %>%
  unnest(data) %>%
  rename(alteration_group=category) %>%
  mutate(confounder="none") %>%
  group_by(cancer_group, alteration_group, confounder) %>%
  nest()

singlecancer <- combined %>%
  ungroup() %>%
  filter(group != "overall") %>%
  unnest(data) %>%
  rename(alteration_group=category) %>%
  mutate(cancer_group="single") %>%
  rename(confounder=group)
tmp <- filter(singlecancer, confounder=="Inactivating TMB") %>%
  mutate(confounder="TMB")
singlecancer <- filter(singlecancer, confounder != "Inactivating TMB") %>%
  bind_rows(tmp)
singlecancer2  <- singlecancer %>%
  group_by(cancer_group, alteration_group, confounder) %>%
  nest()

combined2 <- bind_rows(multicancer, singlecancer2)
saveRDS(combined2,
        file=here(outdir, "combined.rds"))

dat <- combined2 %>%
  filter(cancer_group=="multi") %>%
  unnest("data")
pal <- dat %>%
  filter(p > cutoff_ntest | abs(`50%`) > 6) %>%
  group_by(cancer) %>%
  summarize(n=n(), .groups="drop") %>%
  arrange(-n) %>%
  select(-n)
show_col(viridis_pal()(nrow(pal)))
display.brewer.all(n=nrow(pal), type="all", select=NULL, exact.n=TRUE)
pal$color <- brewer.pal(n=nrow(pal), "Dark2")
saveRDS(pal, file=here(outdir, "colors.rds"))





