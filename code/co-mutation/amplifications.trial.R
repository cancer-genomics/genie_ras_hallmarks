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
outdir <- here("output", "co-mutation", "amplifications.R")
fs::dir_create(outdir)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
library(magrittr)
library(tidyverse)
library(dplyr)
library(contingency.table)
source(here("code", "co-mutation", "stan_funs.R"))
data(amplification_ct, package="contingency.table")

print(paste0("index is ", INDEX) )


mt  <- amplification_ct %>%
    unnest(data)  %>%
    group_by(category, nonras) %>%
    nest()

print("done mt")

gene <- unique(mt$nonras)[INDEX]

print(gene)
mt2 <- mt %>%
    filter(nonras==gene)
print(dim(mt2))
print("done mt2")

stanfile <- here("code", "stan", "log_lin_hierarchical.stan")
print("done stanfile")
model <- stan_model(stanfile)
print("done model") 


#params <- list(chains=3, thin=10, iter=10e3, warmup=5000,
               stanmodel=model)

#multicancer <- fit_stan("overall", mt2, params)
#print("done multicancer")
#sex <- fit_stratify("sex", mt2, params)


