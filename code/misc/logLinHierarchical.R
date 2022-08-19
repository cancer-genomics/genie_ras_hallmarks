library(magrittr)
library(tidyverse)
library(genie.data)
library(rasfuns)
library(fs)
library(rjags)
library(ggmcmc)
library(rstan)
library(bridgesampling)
library(MASS)
library(ggpubr)
library(qdapRegex)
library(tidybayes)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
data("mutation_data-6_1", package="genie.data")

gene <- "K"
codons <- "12"
study <- "matched"
min_N <- 10
analyses <-
  expand.grid(gene, codons, study, min_N) %>%
  as_tibble() %>%
  set_colnames(c("gene", "codons", "study", "min_samples"))


### Get matched data
matched.dat <- mutationDataWithMatchedNormal(mutation_data)

nsclc <- "Non-Small Cell Lung Cancer"


kras_tp53_dat <-
  matched.dat %>%
  mutate(codon_start=codonStart3(.)) %>%
  filter(!is.na(codon_start))

amino_acid <- translateNucleotides(kras_tp53_dat$codons[kras_tp53_dat$hugo_symbol == "KRAS" & kras_tp53_dat$codon_start == 12])
kras_tp53_dat$amino_acid <- ifelse(kras_tp53_dat$hugo_symbol == "KRAS" & kras_tp53_dat$codon_start == 12, amino_acid, NA)

kras12c_tp53_dat <-
  kras_tp53_dat %>%
  group_by(tumor_sample_barcode, cancer_type) %>%
  summarize(kras12 = sum(codon_start == 12 & hugo_symbol == "KRAS" & amino_acid == "G/C") >= 1,
            tp53 = sum(hugo_symbol == "TP53") >= 1,
            egfr = sum(hugo_symbol == "EGFR") >= 1) %>%
  ungroup() %>%
  mutate_all(.funs = as.factor) %>%
  group_by(cancer_type, kras12, tp53, .drop = FALSE) %>%
  tally() %>%
  ungroup() %>%
  rename(y = n) %>%
  group_by(cancer_type) %>%
  mutate(n = sum(y)) %>%
  ungroup() %>%
  arrange(desc(-n)) %>%
  filter(n >= min_N)

kras_tp53_dat <-
  kras_tp53_dat %>%
  group_by(tumor_sample_barcode, cancer_type) %>%
  summarize(kras12 = sum(codon_start == 12 & hugo_symbol == "KRAS") >= 1,
            tp53 = sum(hugo_symbol == "TP53") >= 1,
            egfr = sum(hugo_symbol == "EGFR") >= 1) %>%
  ungroup() %>%
  mutate_all(.funs = as.factor) %>%
  group_by(cancer_type, kras12, tp53, .drop = FALSE) %>%
  tally() %>%
  ungroup() %>%
  rename(y = n) %>%
  group_by(cancer_type) %>%
  mutate(n = sum(y)) %>%
  ungroup() %>%
  arrange(desc(-n)) %>%
  filter(n >= min_N)


##########################################
### hierarchical model
log_lin_hierarch <- stan_model('log_lin_hierarchical.stan')

### Create matrix for each cancer type
cancers <- unique(kras_tp53_dat$cancer_type)
C <- length(cancers)
ngenes <- 2
y <- matrix(NA, nrow = C, ncol = ngenes ^ 2)
### Create model matrix
x_df <- expand.grid(g1 = c(F, T), g2 = c(F, T))
xmat <- model.matrix(~g1*g2, data = x_df)

### Now create y
for(i in 1:nrow(y)) {
  for(j in 1:ncol(y)) {
    cancer_df <- kras_tp53_dat %>% filter(cancer_type == cancers[i])
    y[i,j] <- cancer_df$y[cancer_df$kras12 == x_df$g1[j] & cancer_df$tp53 == x_df$g2[j]]
  }
}

data_saturated <- list(y = y, N = ncol(y), C = C, x = xmat, K = ncol(xmat))
### Run hierarchical model for saturated model
interaction_model <- sampling(log_lin_hierarch, data = data_saturated, iter = 5000)

### Get interaction effects
interaction_betas <-
  ggs(interaction_model, family = "^beta") %>%
  filter(grepl(",4", Parameter)) %>%
  mutate(cancer_index = unlist(rm_between(Parameter, "[", ",", extract = TRUE)),
         cancer_type = cancers[as.numeric(cancer_index)])

### CI width plot
samp_sizes <- kras_tp53_dat %>% group_by(cancer_type) %>% summarize(n = sum(y))
interaction_betas$cancer_type <- factor(interaction_betas$cancer_type,
                                        levels = samp_sizes$cancer_type[order(samp_sizes$n)],
                                        labels = paste(samp_sizes$cancer_type[order(samp_sizes$n)],
                                                                              "n =", samp_sizes$n[order(samp_sizes$n)]))
ci_tp53_kras12 <-
  interaction_betas %>%
  group_by(cancer_type) %>%
  median_qi(value) %>%
  ggplot(aes(x = value, y = cancer_type)) +
  geom_pointintervalh() +
  xlab("Beta interaction") +
  ylab("") +
  theme(axis.text.y = element_blank())

### Posteriors for posterior probabilities
prev_interaction <-
  ggs(interaction_model, family = "pi") %>%
  mutate(cancer_index = unlist(rm_between(Parameter, "[", ",", extract = TRUE)),
         cancer_type = cancers[as.numeric(cancer_index)],
         prev_index = unlist(rm_between(Parameter, ",", "]", extract = TRUE)))
prev_marginals <-
  prev_interaction %>%
  group_by(Iteration, Chain, cancer_type) %>%
  summarize(kras12c = sum(value[prev_index == 2] + value[prev_index == 4]),
            tp53 = sum(value[prev_index == 3] + value[prev_index == 4]))

prev_marginals$cancer_type <- factor(prev_marginals$cancer_type,
                                     levels = samp_sizes$cancer_type[order(samp_sizes$n)],
                                     labels = paste(samp_sizes$cancer_type[order(samp_sizes$n)],
                                                    "n =", samp_sizes$n[order(samp_sizes$n)]))

kras12c_posterior_prevalence <-
  prev_marginals %>%
  group_by(cancer_type) %>%
  median_qi(kras12c) %>%
  ggplot(aes(x = kras12c, y = cancer_type)) +
  geom_pointintervalh() +
  xlab("kras12c posterior prevalence") +
  ylab("Cancer Type")

library(gtable)
g <- ggplotGrob(kras12c_posterior_prevalence)
s <- gtable_filter(g, 'axis-l|ylab', trim=F)  # use trim depending on need
kras12c_posterior_prevalence <-
  kras12c_posterior_prevalence +
  theme(axis.text.y = element_blank()) +
  ylab("")

tp53_posterior_prevalence <-
  prev_marginals %>%
  group_by(cancer_type) %>%
  median_qi(tp53) %>%
  ggplot(aes(x = tp53, y = cancer_type)) +
  geom_pointintervalh() +
  xlab("TP53 posterior prevalence") +
  ylab("") +
  theme(axis.text.y = element_blank())

kras_tp53_combined_posteriors <- ggarrange(s, kras12c_posterior_prevalence,
                                           tp53_posterior_prevalence,
                                           ci_tp53_kras12,
                                           nrow = 1)


if(!dir.exists("../eda_figs")){
  dir.create("../eda_figs", recursive = TRUE)
}

ci_tp53_kras12 <-
  interaction_betas %>%
  group_by(cancer_type) %>%
  median_qi(value) %>%
  ggplot(aes(x = value, y = cancer_type)) +
  geom_pointintervalh() +
  xlab("Beta interaction") +
  ylab("Cancer Type")
ggsave(file.path("../eda_figs", "ci_kras_12c_tp53_by_cancer.pdf"),
       ci_tp53_kras12, width = 10, height = 6)

ggsave(file.path("../eda_figs", "kras_12c_tp53_three_panels.pdf"),
       kras_tp53_combined_posteriors , width = 12, height = 6)

### Plot observed versus posterior
obs_rates <-
  kras_tp53_dat %>%
  mutate(kras12 = as.logical(kras12), tp53 = as.logical(tp53)) %>%
  mutate(obs = y / n,
         prev_index = ifelse(!kras12 & !tp53, 1,
                             ifelse(kras12 & !tp53, 2,
                                    ifelse(!kras12 & tp53, 3, 4))))
prev_interaction_summary <-
  prev_interaction %>%
  mutate(prev_index = as.numeric(prev_index)) %>%
  group_by(cancer_type, prev_index) %>%
  median_qi(value)

obs_vs_est <- inner_join(prev_interaction_summary, obs_rates)

ggplot(obs_vs_est , aes(x = obs, y = value)) +
  geom_point() +
  geom_segment(aes(x = obs, xend = obs, y = .lower, yend = .upper)) +
  geom_abline(alpha = .6, col = 'grey') +
  facet_wrap(~prev_index, scales = 'free')

### Model with no pooling
log_lin_nopool <- stan_model("log_lin_nopool.stan")
interaction_model_nopool <- sampling(log_lin_nopool, data = data_saturated, iter = 5000)

interaction_betas_nopool <-
  ggs(interaction_model_nopool, family = "^beta") %>%
  filter(grepl(",4", Parameter)) %>%
  mutate(cancer_index = unlist(rm_between(Parameter, "[", ",", extract = TRUE)),
         cancer_type = cancers[as.numeric(cancer_index)])
interaction_betas_nopool$cancer_type <- factor(interaction_betas_nopool$cancer_type,
                                        levels = samp_sizes$cancer_type[order(samp_sizes$n)],
                                        labels = paste(samp_sizes$cancer_type[order(samp_sizes$n)],
                                                       "n =", samp_sizes$n[order(samp_sizes$n)]))

### Plot median estimates
interaction_betas_compare <- bind_rows(mutate(interaction_betas, method = 'hierarchical'),
                                       mutate(interaction_betas_nopool, method = 'No pooling')) %>%
  group_by(Parameter, cancer_type, method) %>%
  summarize(median = median(value))
interaction_betas_compare <- pivot_wider(interaction_betas_compare, names_from = method, values_from = median)

ggplot(interaction_betas_compare, aes(x = hierarchical, y = `No pooling`)) +
  geom_point() +
  geom_abline()

### Full posterior distributions
interaction_betas_compare <- bind_rows(mutate(interaction_betas, method = 'hierarchical'),
                                       mutate(interaction_betas_nopool, method = 'No pooling')) %>%
  group_by(Parameter, cancer_type, method) %>%
  median_qi(value)
hierarch_vs_nopool_plot <-
  ggplot(interaction_betas_compare, aes(x = value, y = cancer_type, color = method)) +
  geom_pointintervalh(alpha = .4)
ggsave(file.path("../eda_figs", "kras_12c_tp53_hierarch_vs_nopool.pdf"),
       hierarch_vs_nopool_plot, width = 10, height = 6)


prev_interaction_nopool <-
  ggs(interaction_model_nopool, family = "pi") %>%
  mutate(cancer_index = unlist(rm_between(Parameter, "[", ",", extract = TRUE)),
         cancer_type = cancers[as.numeric(cancer_index)],
         prev_index = unlist(rm_between(Parameter, ",", "]", extract = TRUE)))


### Independence model
data_independent <- list(y = y, N = ncol(y), C = C, x = xmat[,1:3], K = ncol(xmat[,1:3]))
independence_model <- sampling(log_lin_hierarch, data = data_independent, iter = 5000)

### Compare posterior probability for joint occurence


### Overall means
mu_beta <- ggs(log_lin_samples, family = "mu_beta")
ggs_density(mu_beta)
### Sigmas
sigmas <- ggs(log_lin_samples, family = "sigma")
ggs_density(sigmas)


ggplot(interaction_betas, aes(x = Iteration, y = value, color = factor(Chain))) +
  geom_line() +
  facet_wrap(~cancer_type, scales = "free")


interaction_zscore <-
  interaction_betas %>%
  group_by(cancer_type) %>%
  summarize(z_score = mean(value) / sd(value)) %>%
  mutate(direction = ifelse(z_score >= 0, "Positive correlation", "Negative correlation"))

interaction_zscore$cancer_type <-  reorder(interaction_zscore$cancer_type, -abs(interaction_zscore$z_score))

kras_12c_tp53_by_cancer <-
  ggplot(interaction_zscore, aes(x = cancer_type, y = abs(z_score), color = direction)) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90))
if(!dir.exists("../eda_figs")){
  dir.create("../eda_figs", recursive = TRUE)
}
ggsave(file.path("../eda_figs", "kras_12c_tp53_by_cancer.pdf"),
       kras_12c_tp53_by_cancer, width = 10, height = 6)

