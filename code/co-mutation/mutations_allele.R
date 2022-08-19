
## running co-mut model for  mutation data at  allele level
## source script:: code/co-mutation/data-derived/mutations.R

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
outdir <- here("output", "co-mutation", "mutations_allele.R")
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
data(mutation_allele_ct, package="contingency.table")
mt  <- mutation_allele_ct %>%
  unnest(data)  %>%
  group_by(category, nonras) %>%
  nest()
gene <- unique(mt$nonras)[INDEX]
mt2 <- mt %>%
  filter(nonras==gene)
stanfile <- here("code", "stan", "log_lin_hierarchical.stan")
model <- stan_model(stanfile)
params <- list(chains=3, thin=10, iter=10e3, warmup=5000,
               stanmodel=model)
##params <- list(chains=1, thin=1, iter=1000, warmup=500, stanmodel=model)
multicancer <- fit_stan("overall", mt2, params)
sex <- fit_stratify("sex", mt2, params)
race <- fit_stratify("race", mt2, params)
stype <- fit_stratify("sample_type", mt2, params)
age <- fit_stratify("agecat", mt2, params)
mt2$mcmc <- list(multicancer, sex, race, stype, age)
mt3 <- mt2 %>%
  select(-data) %>%
  unnest(mcmc)
saveRDS(mt3, file=here(outdir, paste0(gene, ".rds")))
stop("done")
q('no')

