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
data(known_drivers, package="rasfuns")
outdir <- here("output", "co-mutation", "summarize_models.R")
fs::dir_create(outdir)
## co-mutation tables
nonras_genes <- function(ct){
    ct$data[[1]] %>%
        pull(nonras) %>%
        unique()
}
metadata <- tibble(script=c("mutations.R",
                            "inactivating.R",
                            "deletions.R",
                            "any_deletion.R",
                            "amplifications.R",
                            "any_amplification.R",
                            "mutations_tmb.R",
                            "inactivating_tmb.R",
                            "mutations_mutsig.R",
                            "rearrangements.R"),
                   category=c("Any mutation",
                              "Inactivating mutation",
                              "Deletion",
                              "Any deletion",
                              "Amplification",
                              "Any amplification",
                              "TMB",
                              "Inactivating TMB",
                              "Signature",
                              "Rearrangement"))
rearrangement_ct2  <- rearrangement_ct %>%
    unite(uid, c(ras, nonras, cancer), sep="~", remove=FALSE) %>%
    group_by(uid) %>%
    nest()
metadata$object <- list(mutation_ct,
                        inactivating_ct,
                        deletion_ct,
                        any_deletion_ct,
                        amplification_ct,
                        any_amplification_ct,
                        mutation_tmb,
                        inactivating_tmb,
                        mutation_msig,
                        rearrangement_ct2)
rename <- dplyr::rename
mt <- metadata[1:6, ] %>%
    pmap_dfr(harmonize_stats, drivers=known_drivers) %>%
    unnest(stan)

mt <- metadata[1, ] %>%
    pmap_dfr(harmonize_stats, drivers=known_drivers) %>%
    unnest(stan)
mt2 <- mt %>%
    group_by(group, ras, nonras) %>%
    nest()
mt3 <- metadata[c(7, 8, 9, 10), ] %>%
    pmap_dfr(harmonize_stats2, drivers=known_drivers) %>%
    unnest(stan) %>%
    mutate(category=case_when(category=="TMB" ~ "Any mutation",
                              category=="Signature"~"Any mutation",
                              category=="Inactivating TMB"~"Inactivating mutation",
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
