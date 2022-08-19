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
mt  <- amplification_ct %>%
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
##params <- list(chains=1, thin=1, iter=1000, warmup=500,
##               stanmodel=model)
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



if(!exists("INDEX")){
    if(length(commandArgs(trailingOnly=TRUE)) > 0){
        args <- commandArgs(trailingOnly=TRUE)
        INDEX <- as.numeric(args[[1]])
    } else {
        INDEX <- 1
    }
}
set.seed(INDEX)
library(fs)
outdir <- here("output",
               "co-mutation",
               "amplifications.R")
dir_create(outdir)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
library(MASS)
library(magrittr)
library(tidyverse)
library(genie.6.1)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggmcmc)
library(dplyr)
devtools::load_all(here("code", "rasfuns"))
devtools::load_all(here("code", "contingency.table"))
data(amplifications, package="contingency.table")
ras_genes <- rownames(ras_all)
all_genes <- colnames(amplifications)[24:ncol(amplifications)]
nonras_gene <- all_genes[INDEX]
alterations <- amplifications
##
## List of length 20: 1 contingency table per RAS variant
##
ct <- co_occurrence(mutations=alterations, "RAS_12-61", "FH") %>%
    filter(cancer=="non-small cell lung cancer")


ct.list <- ras_genes %>%
    map(co_occurrence, nonras_gene=nonras_gene,
        mutations=alterations) %>%
    setNames(ras_genes)
stan.data <- ct.list %>%
    map(stanCooccurData)
stan.data2 <- stan.data[ !sapply(stan.data, is.null) ]
params <- list(chains=3, thin=2, iter=5e3, warmup=2500)
stanfile <- file.path("..", "stan", "log_lin_hierarchical.stan")
log_lin_hierarch <- stan_model(stanfile)
sampler <- function(stan.data, model, params, ras_gene,
                    nonras_gene){
    samples <- sampling(model,
                        data=stan.data,
                        iter=params$iter,
                        thin=params$thin,
                        chains=params$chains,
                        warmup=params$warmup)
    probs <- summarize_model(samples,
                             as.character(ras_gene),
                             as.character(nonras_gene),
                             stan.data)
    probs2 <- probs %>%
        select(c(ras_gene, nonras_gene, group, Parameter,
                 colnames(.)))
    probs2    
}
results <- stan.data2 %>%
    map2(.y=names(stan.data2),
         .f=sampler, model=log_lin_hierarch,
         params=params, nonras_gene=nonras_gene)
fname <- file.path(outdir, paste0(nonras_gene, ".rds"))
ct.list2 <- ct.list[ names(results) ]
ct.list3 <- ct.list2 %>%
    map(function(x){
        x2 <- x %>%
            unite("comutation", c(ras, mutation), sep="") %>%
            pivot_wider(names_from="comutation", values_from="n")
        rsums <- rowSums(select(x2, -cancer))
        x2$n <- rsums
        ## use same cutoff as `stanCooccurData`
        ##filter(x2, n >= 20)
        x2
    })
L <- map_int(ct.list3, nrow)
comutation_id <- paste(names(results), nonras_gene, sep="~")
cts <- do.call(bind_rows, ct.list3) %>%
    mutate(uid=rep(comutation_id, L))
results2 <- do.call(bind_rows, results) %>%
    unite("uid", c(ras_gene, nonras_gene), sep="~") %>%
    rename(cancer=group) %>%
    left_join(cts, by=c("uid", "cancer"))
saveRDS(results2, file=fname)    


