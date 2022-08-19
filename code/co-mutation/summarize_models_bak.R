library(devtools)
library(magrittr)
library(tidyverse)
library(broom)
library(contingency.table)
library(here)
source(here("code", "co-mutation", "stan_funs.R"))
data(known_drivers, package="rasfuns")
outdir <- here("output", "summarize-com", "tables.R")
fs::dir_create(outdir)

## prevalence



## co-mutation tables
nonras_genes <- function(ct){
    ct$data[[1]] %>%
        pull(nonras) %>%
        unique()
}
metadata <- tibble(script=c("mutations.R",
                            "inactivating.R",
                            "deletions.R",
                            "amplifications.R",
                            "mutations_tmb.R",
                            "inactivating_tmb.R",
                            "mutations_mutsig.R",
                            "rearrangements.R",
                            "any_deletion.R",
                            "any_amplification.R"),
                   category=c("Any mutation",
                              "Inactivating mutation",
                              "Deletion",
                              "Amplification",
                              "TMB",
                              "Inactivating TMB",
                              "Signature",
                              "Rearrangement",
                              "Any deletion",
                              "Any amplification"))
rearrangement_ct2  <- rearrangement_ct %>%
    unite(uid, c(ras, nonras, cancer), sep="~", remove=FALSE) %>%
    group_by(uid) %>%
    nest()
metadata$object <- list(mutation_ct,
                        inactivating_ct,
                        deletion_ct,
                        amplification_ct,
                        mutation_tmb,
                        inactivating_tmb,
                        mutation_msig,
                        rearrangement_ct2,
                        any_deletion_ct,
                        any_amplification_ct)
mt <- metadata[c(1:4, 9:10), ] %>%
    pmap_dfr(harmonize_stats, drivers=known_drivers) %>%
    unnest(stan)
mt2 <- mt %>%
    group_by(group, ras, nonras) %>%
    nest()
mt3 <- metadata[5:8, ] %>%
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
saveRDS(combined,
        file=here("output", "summarize-com", "tables.R",
                  "results.rds"))

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
        file=here("output", "summarize-com", "tables.R",
                  "combined.rds"))
stop("done")
q('no')



## Inactivating mutations
data(inactivating_ct, package="contingency.table")
nonrasgenes <- nonras_genes(inactivating_ct)
mt2 <- here("output", "co-mutation", "inactivating.R",
            paste0(nonrasgenes, ".rds")) %>%
    map_dfr(readRDS) %>%
    rename(group=category) %>%
    mutate(category="Any mutation") %>%
    group_by(group, ras) %>%
    nest()
mt2$stats <- mt2$data %>%
    map(p_adjust, known_drivers)

inactivating <- parse_stan(mt.list, variant="Inactivating mutation",
                           known_drivers)

data(deletion_ct, package="contingency.table")
rasgenes <- ras_genes(deletion_ct)
mt.list <- here("output", "co-mutation", "deletions.R",
           paste0(rasgenes, ".rds")) %>%
    map(readRDS)
deletion <- parse_stan(mt.list, variant="Deletion",
                       known_drivers)

data(amplification_ct, package="contingency.table")
rasgenes <- ras_genes(amplification_ct)
mt.list <- here("output", "co-mutation", "amplifications.R",
           paste0(rasgenes, ".rds")) %>%
    map(readRDS)
amplification <- parse_stan(mt.list, variant="Amplification",
                            known_drivers)

if(FALSE){
    data(rearrangement_ct, package="contingency.table")
    rasgenes <- unique(rearrangement_ct$ras)
    mt.list <- here("output", "co-mutation", "rearrangements.R",
                    paste0(rasgenes, ".rds")) %>%
        map(readRDS)
    rearrangement <- parse_stan(mt.list, variant="Rearrangement",
                                known_drivers)

    data(mutation_tmb, package="contingency.table")
    mutation_tmb2 <- mutation_tmb %>%
        separate(uid, c("ras", "nonras", "cancer"), sep="~")
    genes <- unique(mutation_tmb2$nonras)
    mt.list <- here("output", "co-mutation", "mutations_tmb.R",
               paste0(genes, ".rds")) %>%
        map(readRDS)
    parse_tmb <- function(x, drivers){
        mt <- x
        mt2 <- mt$mcmc[[1]]
        mt3 <- unnest(mt2, mcmc)
        dat <- unnest(mt2, data) %>%
            mutate(Parameter=factor(variable),
                   Parameter=as.integer(Parameter),
                   Parameter=paste0("beta[", Parameter, ",4]"))
        mt4 <- mt3 %>%
            select(-data) %>%
            left_join(select(dat, -variable), by=c("cancer", "Parameter")) %>%
            ungroup()
        results <- p_adjust(mt4, drivers=drivers)
    }
    mt_tmb <- mt.list %>%
        map_dfr(parse_tmb, known_drivers) %>%
        mutate(category="TMB", variant="Any mutation") %>%
        group_by(category, variant) %>%
        nest()

    data(inactivating_tmb, package="contingency.table")
    mutation_tmb2 <- inactivating_tmb %>%
        separate(uid, c("ras", "nonras", "cancer"), sep="~")
    genes <- unique(mutation_tmb2$ras)
    mt.list <- here("output", "co-mutation", "inactivating_tmb.R",
                    paste0(genes, ".rds")) %>%
        map(readRDS)
    imt_tmb <- mt.list %>%
        map_dfr(parse_tmb, known_drivers) %>%
        mutate(category="TMB", variant="Inactivating mutation") %>%
        group_by(category, variant) %>%
        nest()
}
mt <- list(mutation,
     inactivating,
     deletion,
     amplification) %>%
    map_dfr(bind_rows)
tmp <- mt$data %>%
    map(function(x) select(x, -category))
mt$data <- tmp
mt2 <- unnest(mt, data)
saveRDS(mt2, here("output", "summarize-com", "tables.R",
                  "results.rds"))

co_mutation <- mtable %>%
    bind_rows(itable) %>%
    bind_rows(dtable) %>%
    bind_rows(atable) %>%
    separate(uid, c("ras", "nonras"), sep="~",
             remove=FALSE) %>%
    group_by(ras, type) %>%
    nest()
saveRDS(co_mutation, file=file.path(outdir, "co_mutation.rds"))

