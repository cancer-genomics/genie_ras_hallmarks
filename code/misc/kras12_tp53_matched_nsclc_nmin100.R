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
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
data("mutation_data-6_1", package="genie.data")

gene <- c("K")
codons <- c("12")
study <- c("matched")
min_N <- 100
gene <- "K"
codons <- "12"
study <- "matched"
min_N <- 100
analyses <- expand.grid(gene, codons, study, min_N) %>%
  as_tibble() %>%
  set_colnames(c("gene", "codons", "study", "min_samples"))


### Get matched data
matched.dat <- mutationDataWithMatchedNormal(mutation_data)

nsclc <- "Non-Small Cell Lung Cancer"


#### Get data for all cancers
kras_tp53_dat <-
  matched.dat %>%
  mutate(codon_start=codonStart3(.)) %>%
  filter(!is.na(codon_start)) %>%
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
  arrange(desc(-n))


###### Log-linear model
nsclc_dat <-
  kras_tp53_dat %>%
  filter(cancer_type == nsclc)
x_ind <- model.matrix(y ~ kras12 + tp53, data = nsclc_dat)
x_sat <- model.matrix(y ~ kras12*tp53, data = nsclc_dat)

log_lin_model <- stan_model('log_lin.stan')
log_lin_indep <- sampling(log_lin_model, data = list(y = nsclc_dat$y,
                                                 x = x_ind,
                                                 N = nrow(nsclc_dat),
                                                 K = ncol(x_ind)))
log_lin_dep <- sampling(log_lin_model, data = list(y = nsclc_dat$y,
                                                    x = x_sat,
                                                    N = nrow(nsclc_dat),
                                                    K = ncol(x_sat)))
beta_indep <- ggs(log_lin_indep, family = "beta")
beta_dep <- ggs(log_lin_dep, family = "beta")
ggs_density(beta_indep)
ggs_density(beta_dep)
pi_indep <- ggs(log_lin_indep, family = "pi")
pi_dep <- ggs(log_lin_dep, family = "pi")
ggs_density(pi_indep)
ggs_density(pi_dep)



### Try it for cancer with smaller sample size
sellar <- subset(kras_tp53_dat, cancer_type == "Sellar Tumor")
x_sat <- model.matrix(y ~ kras12*tp53, data = sellar)

log_lin_dep <- sampling(log_lin_model, data = list(y = sellar$y,
                                                   x = x_sat,
                                                   N = nrow(sellar),
                                                   K = ncol(x_sat)))
sigma2_sellar <- ggs(log_lin_dep, family = "sigma2")
beta_sellar <- ggs(log_lin_dep, family = "beta")

##########################################
### hierarchical model
log_lin_hierarch <- stan_model('log_lin_hierarchical.stan')
### Get data for two biggest cancer types
breast_nsclc <- filter(kras_tp53_dat, cancer_type %in% c(nsclc, "Breast Cancer"))
cancers <- unique(breast_nsclc$cancer_type)
C <- length(cancers)
y_all <- matrix(NA, nrow = C, ncol = 4)
for(c in 1:C) {
  y_all[c,] <- breast_nsclc %>% filter(cancer_type == cancers[c]) %>% pull(y)
}

hierarch_samples <- sampling(log_lin_hierarch, data = list(y = y_all, x = x_sat,
                             N = ncol(y_all), K = nrow(x_sat), C = C))
betas <- ggs(hierarch_samples, family = "^beta")
ggs_density(betas) +
  facet_wrap(~Parameter, nrow = 2, scales = "free")

pis <- ggs(hierarch_samples, family = "pi")
ggs_density(pis) +
  facet_wrap(~Parameter, nrow = 2, scales = "free")


############# Now add in leukemia
breast_nsclc_leuk <- filter(kras_tp53_dat,
                            cancer_type %in% c(nsclc, "Breast Cancer", "Leukemia"))
cancers <- unique(breast_nsclc_leuk$cancer_type)
C <- length(cancers)
y_all <- matrix(NA, nrow = C, ncol = 4)
for(c in 1:C) {
  y_all[c,] <- breast_nsclc_leuk %>% filter(cancer_type == cancers[c]) %>% pull(y)
}

hierarch_samples <- sampling(log_lin_hierarch,
                             data = list(y = y_all, x = x_sat,
                                         N = ncol(y_all), K = nrow(x_sat), C = C))
betas <- ggs(hierarch_samples, family = "^beta")
ggs_density(betas) +
  facet_wrap(~Parameter, nrow = C, scales = "free")

pis <- ggs(hierarch_samples, family = "pi")
ggs_density(pis) +
  facet_wrap(~Parameter, nrow = C, scales = "free")

########################################
### Add in hodgkin
breast_nsclc_leuk_hodg <- filter(kras_tp53_dat,
                            cancer_type %in% c(nsclc, "Breast Cancer", "Leukemia",
                                               "Hodgkin Lymphoma"))
cancers <- unique(breast_nsclc_leuk_hodg$cancer_type)
C <- length(cancers)
y_all <- matrix(NA, nrow = C, ncol = 4)
for(c in 1:C) {
  y_all[c,] <- breast_nsclc_leuk_hodg %>% filter(cancer_type == cancers[c]) %>% pull(y)
}

hierarch_samples <- sampling(log_lin_hierarch,
                             data = list(y = y_all, x = x_sat,
                                         N = ncol(y_all), K = nrow(x_sat), C = C))
mu_beta <- ggs(hierarch_samples, family = "mu_beta")
ggs_density(mu_beta) +
  facet_wrap(~Parameter, nrow = C, scales = "free")

betas <- ggs(hierarch_samples, family = "^beta")
ggs_density(betas) +
  facet_wrap(~Parameter, nrow = C, scales = "free")

pis <- ggs(hierarch_samples, family = "pi")
ggs_density(pis) +
  facet_wrap(~Parameter, nrow = C, scales = "free")


########################################
### KRAS-TP53-EGFR in NSCLC
kras_tp53_egfr <-
  matched.dat %>%
  mutate(codon_start=codonStart3(.)) %>%
  filter(!is.na(codon_start)) %>%
  group_by(tumor_sample_barcode, cancer_type) %>%
  summarize(kras12 = sum(codon_start == 12 & hugo_symbol == "KRAS") >= 1,
            tp53 = sum(hugo_symbol == "TP53") >= 1,
            egfr = sum(hugo_symbol == "EGFR") >= 1) %>%
  ungroup() %>%
  mutate_all(.funs = as.factor) %>%
  group_by(cancer_type, kras12, tp53, egfr, .drop = FALSE) %>%
  tally() %>%
  ungroup() %>%
  rename(y = n) %>%
  group_by(cancer_type) %>%
  mutate(n = sum(y)) %>%
  ungroup() %>%
  arrange(desc(-n))

kras_tp53_egfr_nsclc <- filter(kras_tp53_egfr, cancer_type == nsclc)
x_sat <- model.matrix(y ~ kras12*tp53*egfr, data = kras_tp53_egfr_nsclc)
kras_tp53_egfr_nsclc_br <- filter(kras_tp53_egfr, cancer_type %in% c(nsclc, "Breast Cancer"))
cancers <- unique(kras_tp53_egfr_nsclc_br$cancer_type)
C <- length(cancers)
y_all <- matrix(NA, nrow = C, ncol = 8)
for(c in 1:C) {
  y_all[c,] <- kras_tp53_egfr_nsclc_br %>% filter(cancer_type == cancers[c]) %>% pull(y)
}
hierarch_samples <- sampling(log_lin_hierarch,
                             data = list(y = y_all, x = x_sat,
                                         N = ncol(y_all), K = nrow(x_sat), C = C))
betas <- ggs(hierarch_samples, family = "^beta")
ggs_density(betas) +
  facet_wrap(~Parameter, nrow = C, scales = "free")


#######################################

y <- as.matrix(kras_tp53_dat[,c("kras12", "tp53")])

data <- list(Y = structure(as.integer(y), .Dim = c(nrow(y), 2)), N = nrow(y))
stan_probit <- stan(file = 'multivariate_probit.stan', data = data)
rho_samples <- ggs(stan_probit, family = "rho")
rho_plot <- ggs_density(rho_samples)

post_samples <- ggs(stan_probit) %>% filter(!(grepl("log_lik", Parameter)))
post_samples_wide <- pivot_wider(post_samples, names_from = Parameter,values_from = value)

### Note that some have greater than one mutation for each
### Create matrix
### rows will be kras mutation/no-mutation, column will be tp53
comutation_mat <- matrix(NA, nrow = 2, ncol = 2)
rownames(comutation_mat) <- c("kras", "no_kras")
colnames(comutation_mat) <- c('tp53', 'no_tp53')
comutation_mat[1,1] <- sum(kras_tp53_dat$kras12 >= 1 & kras_tp53_dat$tp53 >= 1)
comutation_mat[1,2] <- sum(kras_tp53_dat$kras12 >= 1 & kras_tp53_dat$tp53 == 0)
comutation_mat[2,1] <- sum(kras_tp53_dat$kras12 == 0 & kras_tp53_dat$tp53 >= 1)
comutation_mat[2,2] <- sum(kras_tp53_dat$kras12 == 0 & kras_tp53_dat$tp53 == 0)

post_pred <- lapply(1:nrow(post_samples_wide), function(i) {
  mu <- c(post_samples_wide$`mu[1]`[i], post_samples_wide$`mu[2]`[i])
  Sigma <- diag(1, 2)
  Sigma[1,2] <-Sigma[2,1]<- post_samples_wide$rho[i]
  z <- mvrnorm(n = nrow(y), mu = mu, Sigma = Sigma)
  y_samp <- z >= 0
  tabs <- table(y_samp[,1], y_samp[,2])
  df <- data.frame(mut = c("-/-", "+/-", "-/+", "+/+"),
                   post = as.vector(tabs),
                   truth = c(comutation_mat[2,2], comutation_mat[1,2],
                             comutation_mat[2,1], comutation_mat[1,1]),
                   sample = i)
  return(df)
})
post_pred_df <- do.call(rbind, post_pred)
dependence_plot <-
  ggplot(post_pred_df, aes(x = post)) +
  geom_line(stat = 'density') +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~mut, nrow =2) +
  xlab("Posterior predictive count") +
  ggtitle("Dependence model")

### Now do this for independence
stan_probit_ind <- stan(file = 'multivariate_probit_indep.stan', data = data)

post_samples <- ggs(stan_probit_ind)%>% filter(!(grepl("log_lik", Parameter)))
post_samples_wide <- pivot_wider(post_samples, names_from = Parameter,values_from = value)

post_pred <- lapply(1:nrow(post_samples_wide), function(i) {
  mu <- c(post_samples_wide$`mu[1]`[i], post_samples_wide$`mu[2]`[i])
  Sigma <- diag(1, 2)
  z <- mvrnorm(n = nrow(y), mu = mu, Sigma = Sigma)
  y_samp <- z >= 0
  tabs <- table(y_samp[,1], y_samp[,2])
  df <- data.frame(mut = c("-/-", "+/-", "-/+", "+/+"),
                   post = as.vector(tabs),
                   truth = c(comutation_mat[2,2], comutation_mat[1,2],
                             comutation_mat[2,1], comutation_mat[1,1]),
                   sample = i)
  return(df)
})
post_pred_df <- do.call(rbind, post_pred)
ind_plot <-
  ggplot(post_pred_df, aes(x = post)) +
  geom_line(stat = 'density') +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~mut, nrow =2) +
  xlab("Posterior predictive count") +
  ggtitle("Independence model")

combined_plot <- ggarrange(ggarrange(dependence_plot, ind_plot, ncol = 2),
                           rho_plot, nrow = 2)
### Compare LOO
library("loo")

# Extract pointwise log-likelihood and compute LOO
log_lik_1 <- extract_log_lik(stan_probit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1))
loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)

log_lik_2 <- extract_log_lik(stan_probit_ind, merge_chains = FALSE)
r_eff_2 <- relative_eff(exp(log_lik_2))
loo_2 <- loo(log_lik_2, r_eff = r_eff_2, cores = 2)

comp <- compare(loo_1, loo_2)
print(comp)


### Create vector for comutations (goes by column)
# kras+tp53+, kras-tp53+, kras+tp53-, kras-tp53-
y_dep <- as.vector(comutation_mat)
### Independence vector (kras sum and tp53 sum)
y_indep <- c(sum(comutation_mat[1,]), sum(comutation_mat[,2]))
n <- sum(comutation_mat)

comutation_dep <- stan(file = 'comutation_dep.stan', data = list(y = y_dep))
comutation_indep <- stan(file = 'comutation_indep.stan', data = list(y = y_indep, n = n))

### Look at difference in conditional probabilities
cond_prob <- ggs(comutation_dep, family = "cond_probs")
ggs_caterpillar(cond_prob)

### Compute bayes factor
BayesFactor::contingencyTableBF(comutation_mat, sampleType = "jointMulti")

bayes_factor(bridge_sampler(comutation_dep), bridge_sampler(comutation_indep))
theta_dep <- ggs(comutation_dep)
theta_indep <- ggs(comutation_indep)


