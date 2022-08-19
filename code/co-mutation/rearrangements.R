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
outdir <- here("output", "co-mutation", "rearrangements.R")
fs::dir_create(outdir)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
library(magrittr)
library(tidyverse)
library(genie.6.1)
library(dplyr)
library(contingency.table)
source(here("code", "co-mutation", "stan_funs.R"))
data(rearrangement_ct, package="contingency.table")
mt  <- rearrangement_ct %>%
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
##trace(fit_stratify, browser)
stype <- fit_stratify("sample_type", mt2, params)
##age <- fit_stratify("agecat", mt2, params)
mt2$mcmc <- list(multicancer, sex, race, stype)
mt3 <- mt2 %>%
    select(-data) %>%
    unnest(mcmc)
saveRDS(mt3, file=here(outdir, paste0(gene, ".rds")))
stop("done")
q('no')


mt  <- rearrangement_ct %>%
    ##unnest(data)  %>%
    group_by(category, ras) %>%
    nest()
ras_gene <- unique(mt$ras)[INDEX]
mt2 <- mt %>%
    filter(ras==ras_gene)
stan_format <- function(x){
    y <- select(x, c("00", "10", "01", "11")) %>%
        as.matrix()
    rownames(y) <- x$cancer
    keep <- rowSums(y) >= 20
    y <- y[keep, , drop = FALSE]
    C <- nrow(y)
    N <- ncol(y)
    x <- cbind(1, c(0, 0, 1, 1), c(0, 1, 0, 1), c(0, 0, 0, 1))
    dat.list <- list(y = y, C = C, N = N, x = x, K = ncol(x))
    dat.list
}
sampler <- function(stan.data, params){
    samples <- sampling(params$stanmodel,
                        data=stan.data,
                        iter=params$iter,
                        thin=params$thin,
                        chains=params$chains,
                        warmup=params$warmup)

    nms <- names(rstan::extract(samples))
    pars <- c("marginal",
              "beta",
              "mu_beta[3]",
              "sigma")
    xx <- rstan::summary(samples,
                         pars=pars,
                         probs=c(0.025, 0.05, 0.25, 0.5,
                                 0.75, 0.95, 0.975))$summary
    stats <- as_tibble(xx) %>%
        mutate(Parameter=rownames(xx))
    y <- stan.data$y
    stats2 <- tibble(index=seq_len(nrow(y)),
                      variable=rownames(y)) %>%
        mutate(Parameter=paste0("beta[", index, ",4]")) %>%
        bind_rows(tibble(variable="overall", Parameter="mu_beta[3]")) %>%
        select(-index) %>%
        left_join(stats, by="Parameter")
    stats2
}
fit_stan <- function(x, params) {
    x <- filter(x, number_ras > 25, number_nonras > 25)
    nonras <- x %>%
        group_by(nonras) %>%
        nest()
    stan.data <- nonras$data %>%
        map(stan_format) %>%
        setNames(nonras$nonras)
    nonras$mcmc <- stan.data %>%
        map(sampler,
            params=params)
    nonras
}
stanfile <- here("code", "stan", "log_lin_hierarchical.stan")
model <- stan_model(stanfile)
params <- list(chains=3, thin=10, iter=10e3, warmup=5000,
               stanmodel=model)
mt2$mcmc <- mt2$data %>%
    map(fit_stan, params=params)
saveRDS(mt2, file=here(outdir, paste0(ras_gene, ".rds")))
q('no')

##here::i_am("code/co-mutation/rearrangements.R")
##library(here)
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
outdir <- file.path("..", "..", "output",
                    "co-mutation", "rearrangements.R")
dir_create(outdir)
library(rstan)
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())
library(MASS)
library(magrittr)
library(tidyverse)
library(genie.6.1)
##library(rasfuns)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggmcmc)
library(dplyr)
devtools::load_all("../../../rasfuns")
devtools::load_all("../../../comutation_wflow/code/contingency.table/")
data(rearrangements, package="contingency.table")
ras_gene <- rownames(ras_all)[INDEX]
all_genes <- colnames(rearrangements)[26:ncol(rearrangements)]
tmp <- select(rearrangements, all_of(all_genes))
freq <- colSums(tmp)
all_genes2 <- colnames(tmp)[freq >= 50]
nms <- c(colnames(rearrangements)[1:23],
         all_genes2)
rearrangements2 <- select(rearrangements, all_of(nms))
ct.list <- all_genes2 %>%
    map(co_occurrence, ras_gene=ras_gene,
        mutations=rearrangements) %>%
    setNames(all_genes2)
stan.data <- ct.list %>%
    map(stanCooccurData)
ct.list2 <- map(stan.data, function(x) x$y)
ct.list3 <- ct.list2[ !map_lgl(ct.list2, is.null) ]
## only consider cancers with at least 1 fusion
filter_cancers <- function(x){
    sel <- x[, "01"] > 5 | x[, "11"] > 5
    x[sel, , drop=FALSE]
}
ct.list4 <- ct.list3 %>%
    map(filter_cancers)
ct.list4

## Fisher exact test
fisher_test <- function(x){
    broom::tidy(fisher.test(matrix(x, 2, 2))) %>%
        select(estimate, p.value)
}
fisher_apply <- function(x){
    xx <- apply(x, 1, fisher_test) %>%
        map_dfr(bind_rows) %>%
        mutate(cancer=rownames(x))
}
fe <- ct.list4 %>%
    map_dfr(fisher_apply)
freq <- do.call(rbind, ct.list4)
freq2 <- freq %>%
    as_tibble() %>%
    mutate(cancer=rownames(freq),
           rear_id=rep(names(ct.list4), map_int(ct.list4, nrow)))
freq3 <- bind_cols(freq2, select(fe, -cancer))
freq3
