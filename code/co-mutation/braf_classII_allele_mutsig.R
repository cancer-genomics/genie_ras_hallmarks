
if(!exists("INDEX")){
  if(length(commandArgs(trailingOnly=TRUE)) > 0){
    args <- commandArgs(trailingOnly=TRUE)
    INDEX <- as.numeric(args[[1]])
  } else {
    INDEX <- 1
  }
}
set.seed(INDEX)
library(here)
outdir <- here("output", "co-mutation", "braf_classII_allele_mutsig.R")
fs::dir_create(outdir)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
library(magrittr)
library(tidyverse)
library(genie.6.1)
library(dplyr)
library(purrr)
library(contingency.table)
source(here("code", "co-mutation", "stan_funs.R"))
data(braf_classII_allele_msig, package="contingency.table")
mt  <- braf_classII_allele_msig %>%
  unnest(data)  %>%
  group_by(category, nonras) %>%
  nest()
gene <- unique(mt$nonras)[INDEX]
mt2 <- mt %>%
  filter(nonras==gene)
stanfile <- here("code", "stan", "log_lin_hierarchical.stan")
model <- stan_model(stanfile)
#params <- list(chains=3, thin=10, iter=10e3, warmup=5000,
 #              stanmodel=model)
params <- list(chains=1, thin=1, iter=1000, warmup=500,
               stanmodel=model)
msig <- fit_stratify("signature", mt2, params) %>%
  mutate(nonras=mt2$nonras, category="Mutation signature")
saveRDS(msig, file=here(outdir, paste0(gene, ".rds")))
stop("done")
q('no')



