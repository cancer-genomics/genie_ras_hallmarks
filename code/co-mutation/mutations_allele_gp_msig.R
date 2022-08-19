
## running co-mut model for  mutation data at  allele level
## source script:: code/co-mutation/mutations_msig.R

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
outdir <- here("output", "co-mutation", "mutations_allele_gp_mutsig.R")
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
data(mutation_allele_gp_msig, package="contingency.table")
mt  <- mutation_allele_gp_msig %>%
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
msig <- fit_stratify("signature", mt2, params) %>%
  mutate(nonras=mt2$nonras, category="Mutation signature")
saveRDS(msig, file=here(outdir, paste0(gene, ".rds")))
stop("done")
q('no')


###############################################################
library(rstan)
library(tidyverse)
library(dplyr)
library(fs)
library(magrittr)
library(genie.6.1)
library(purrr)
set.seed(1498)
outputdir <- file.path("..", "..", "output", "tmb")
outdir <- file.path(outputdir, "mutations.R")
dir_create(outdir)
toplevel <- "../../.."
devtools::load_all(file.path(toplevel, "rasfuns"))
devtools::load_all(file.path(toplevel,
                             "comutation_wflow/code",
                             "contingency.table"))
multicancer.path <- file.path("../..",
                              "output", "co-mutation",
                              "mutations.R")
fnames <- file.path(multicancer.path,
                    paste0(rownames(ras_all), ".rds"))
multicancer <- fnames %>%
  "["(file.exists(.)) %>%
  map_dfr(readRDS) %>%
  mutate(z=mean/sd,
         p=ifelse(z < 0, 2*pnorm(z), 2*(1-pnorm(z))),
         p.adj=p.adjust(p, method="BH"))
saveRDS(multicancer,
        file.path(multicancer.path,
                  "multicancer_overall.rds"))
signif <- multicancer %>%
  filter(p.adj < 0.05) %>%
  unite("uid2", c(uid, cancer), sep="_",
        remove=FALSE)
##
##  Add TMB
##
data(mutations, package="contingency.table")
extdir <- system.file("extdata", package="genie.6.1")
fname <- file.path(extdir, "tumor_normal_predictions.csv")
tmb <- read_csv(fname) %>%
  mutate(Estimate=ifelse(Estimate==">64", 64, Estimate)) %>%
  mutate(Estimate=as.numeric(Estimate)) %>%
  mutate(tmb_group=ntile(Estimate, 5),
         tmb_group=ifelse(Estimate >= 40, 6L, tmb_group)) %>%
  select(sample_id, tmb_group) %>%
  filter(sample_id %in% mutations$sample_id)
mutations_tmb <- left_join(tmb, mutations, by="sample_id") %>%
  mutate(tmb_group=factor(tmb_group)) %>%
  filter(tmb_group != 1) %>%
  mutate(tmb_group=droplevels(tmb_group))

## Create contingency table for each of the statistically significant
## associations
signif2 <- signif %>%
  group_by(uid2) %>%
  nest()
ct <- signif2$data %>%
  map_dfr(createCt, mutations_tmb=mutations_tmb) %>%
  group_by(cancer, ras_gene, nonras_gene) %>%
  nest()    

##
## Fit hierarchical log-linear model (model is hierarchical over the TMB groups) 
##
ct$stan.data <- list(ct$data, ct$cancer, ct$ras_gene,
                     ct$nonras_gene) %>%
  pmap(stanFormat)
ct$ras_gene <- as.list(ct$ras_gene)
ct$nonras_gene <- as.list(ct$nonras_gene)

##model <- file.path("..", "stan", "log_lin_hierarchical.stan") %>%
##    stan_model()
model <- file.path("..", "stan", "tmb_model.stan") %>%
  stan_model()
params <- list(chains=3, thin=10, iter=10e3, warmup=5000)
ct$cancer <- as.list(ct$cancer)
results <- ct %>%
  pmap_dfr(sampler, model=model)
ct2 <- ct %>%
  mutate(cancer=unlist(cancer),
         ras_gene=unlist(ras_gene),
         nonras_gene=unlist(nonras_gene)) %>%
  unite("uid2", c("ras_gene", "nonras_gene"), sep="~") %>%
  unite("uid", c(uid2, cancer), sep="_")
ct2$data <- ct2$data %>%
  map(wide_format) 
ct2 <- select(ct2, uid, data) %>%
  set_colnames(c("uid", "counts"))
results2 <- results %>%
  unite("uid2", c("ras_gene", "nonras_gene"), sep="~") %>%
  unite(uid, c("uid2", "cancer"), sep="_") %>%
  group_by(`uid`) %>%
  nest()
results3 <- left_join(results2, ct2, by="uid")
results3$data <- results3$data %>%
  map2(results3$counts, function(x, y){
    tmp <- left_join(x, y, by="group")
  })
results3 <- select(results3, uid, data)
fname <- file.path(outdir, "mutations.rds")
saveRDS(results3, file=fname)

