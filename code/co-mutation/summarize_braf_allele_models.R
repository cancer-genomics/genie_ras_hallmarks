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

outdir <- here("output", "co-mutation", "summarize_braf_allele_models.R")
fs::dir_create(outdir)
## co-mutation tables
nonras_genes <- function(ct){
  ct$data[[1]] %>%
    pull(nonras) %>%
    unique()
}

metadata <- tibble(script=c("braf_classI_allele.R",
                           "braf_classII_allele.R",
			"braf_classIII_allele.R",
			"braf_classI_allele_mutsig.R",
			"braf_classII_allele_mutsig.R",
			"braf_classIII_allele_mutsig.R",
			"braf_classI_allele_tmb.R",
			"braf_classII_allele_tmb.R",
			"braf_classIII_allele_tmb.R"),

                   category=c("BRAF Class I Any mutation",
			      "BRAF Class II Any mutation",
			      "BRAF Class III Any mutation",
			      "BRAF Class I Signature",
                              "BRAF Class II Signature",
                              "BRAF Class III Signature",
			     "BRAF Class I TMB",
                              "BRAF Class II TMB",
                              "BRAF Class III TMB" ))

metadata$object <- list(braf_classI_allele_ct,
                        braf_classII_allele_ct,
			braf_classIII_allele_ct,
			braf_classI_allele_msig,
                        braf_classII_allele_msig,
                        braf_classIII_allele_msig,
			braf_classI_allele_tmb,
                        braf_classII_allele_tmb,
                        braf_classIII_allele_tmb)

rename <- dplyr::rename

mt <- metadata[c(1:3), ] %>%
  pmap_dfr(harmonize_stats, drivers=known_drivers) %>%
  unnest(stan)

mt2 <- mt %>%
  group_by(group, ras, nonras) %>%
  nest()

mt3 <- metadata[-c(1:3), ] %>%
  pmap_dfr(harmonize_stats2, drivers=known_drivers) %>%
 unnest(stan) %>% mutate(category=case_when( category==  "BRAF Class I Signature" ~  "BRAF Class I Any mutation",
                            category==  "BRAF Class II Signature" ~ "BRAF Class II Any mutation",
                             category== "BRAF Class III Signature"~ "BRAF Class III Any mutation",
                            category== "BRAF Class I TMB" ~ "BRAF Class I Any mutation",
                             category== "BRAF Class II TMB" ~ "BRAF Class II Any mutation",
                            category==  "BRAF Class III TMB" ~ "BRAF Class III Any mutation",
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





