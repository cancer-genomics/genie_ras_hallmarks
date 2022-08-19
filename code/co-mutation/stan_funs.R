stan_format <- function(x, var){
    y <- select(x, c("00", "10", "01", "11")) %>%
        as.matrix()
    nms <- select(x, all_of(var))[[1]]
    rownames(y) <- nms
    ##keep <- rowSums(y) >= 20
    ##y <- y[keep, , drop = FALSE]
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

cancer_strata <- function(mc, x){
    overall <- mc %>% map_dfr(function(x){
        filter(x, variable=="overall") %>%
            select(-variable)
    }) %>%
        mutate(ras=x$ras)
    strata <- mc %>% map(function(x) {
        filter(x, variable!="overall") %>%
            select(-variable)
    })
    x$strata <- strata
    strata <- unnest(x, c(data, strata)) %>%
        ungroup()
    total_ <- x$data %>%
        map_dfr(totals)
    total2  <- total_ %>%
        mutate(ras=x$ras)
    overall2 <- left_join(total2, overall, by="ras") %>%
        mutate(variable="overall",
               cancer="multicancer") %>%
        select(colnames(strata))
    results <- bind_rows(strata, overall2)
    results
}

confounder_strata <- function(mc, x){
    overall <- mc %>% map_dfr(function(x){
        filter(x, variable=="overall") %>%
            select(-variable)
    }) %>%
        mutate(cancer=x$cancer,
               ras=x$ras)
    x$strata <- mc %>% map(function(x) {
        filter(x, variable!="overall") %>%
            select(-variable)
    })
    strata <- unnest(x, c(data, strata)) %>%
        ungroup()
    total_ <- x$data %>%
        map_dfr(totals)
    total2  <- total_ %>%
        mutate(cancer=x$cancer,
               ras=x$ras)
    if("uid" %in% colnames(strata)){
        strata <- select(strata, -uid) %>%
            mutate(variable=as.character(variable))
    }
    overall2 <- left_join(total2, overall,
                          by=c("cancer", "ras")) %>%
        mutate(variable="overall") %>%
        select(colnames(strata))
    results <- bind_rows(strata, overall2)
    results
}

totals <- function(x){
    tibble("00" = sum(x$`00`),
           "10" = sum(x$`10`),
           "01" = sum(x$`01`),
           "11" = sum(x$`11`),
           n=sum(x$`n`),
           number_ras=sum(x$number_ras),
           number_nonras=sum(x$number_nonras))
}

fit_stan <- function(category_, x, params) {
    x <- filter(x, category==category_)
    mt3 <- x$data[[1]] %>%
        filter(n >= 20, number_nonras > 10)
    mt4 <- mt3 %>%
        select(-variable) %>%
        group_by(ras) %>%
        nest()
    tmp <- mt4$data %>%
        map(stan_format, "cancer") %>%
        map(sampler, params=params)
    results <- cancer_strata(tmp, mt4)
    results
}

## overall (multiple cancers)
fit_stratify <- function(category_, x, params,
                         MIN_RAS=20, MIN_NONRAS=10) {
    x <- filter(x, category==category_)
    mt3 <- x$data[[1]]  %>%
        filter(n >= MIN_RAS,
               number_nonras > MIN_NONRAS)
    if(nrow(mt3)==0) return(NULL)
    mt4 <- mt3 %>%
        group_by(cancer, ras) %>%
        nest()
    tmp <- mt4$data %>%
        map(stan_format, "variable") %>%
        map(sampler, params=params)
    results <- confounder_strata(tmp, mt4)
    results
}

p_adjust <- function(x, drivers, w=c(0.1, 0.9), alpha=0.05,
                     exclude=c("mu_beta[3]", "sigma"),
                     min_overall=5){ ## minimum marginal frequency
    MIN <- min_overall
    xx <- x %>%
        filter(!Parameter %in% exclude) %>%
        mutate(z=mean/sd,
               p=ifelse(z < 0, 2*pnorm(z), 2*(1-pnorm(z))),
               p=ifelse(p==0, min(p[p > 0]), p)) %>%
        filter(number_ras > MIN, number_nonras > MIN)
    if(nrow(xx) == 0) return(NULL)
    ct <- select(xx, c("00", "10", "01", "11", "n"))
    ct[is.na(ct)] <- 0
    ct$n <- ct$`00` + ct$`10` + ct$`01` + ct$`11`
    chisqfun <- function(x){
        xx <- matrix(x, ncol=2)
        stats <- tidy(chisq.test(xx))
    }
    tmp <- select(ct, -n)
    chisqstats <- apply(tmp, 1, chisqfun) %>%
        map_dfr(bind_rows) %>%
        select(statistic, p.value) %>%
        set_colnames(c("chisq_stat", "chisq_p"))
    ct2 <- bind_cols(ct, chisqstats)
    x2 <- select(xx, -c("00", "10", "01", "11", "n")) %>%
        bind_cols(ct2) %>%
        mutate(chisq_padj=p.adjust(chisq_p, method="BH"))
    ##
    ## Number of unique uids corresponds to number of models evaluated
    ##
    correct_ntests <- x2 %>%
        mutate(is_driver=nonras %in% drivers) %>%
        group_by(is_driver) %>%
        count() %>%
        mutate(labels=as.integer(is_driver) + 1,
               labels=c("non-driver", "driver")[labels]) %>%
        mutate(weight=ifelse(is_driver, w[2], w[1]),
               weight=weight/(sum(weight)),
               alpha=alpha,
               alpha2=alpha * weight/n,
               cutoff_ntest=alpha2/n) %>%
        select(is_driver, cutoff_ntest)
    cutoffs <- x2 %>%
        group_by(nonras) %>%
        summarize(is_driver=nonras %in% drivers,
                  .groups="drop") %>%
        distinct() %>%
        group_by(is_driver) %>%
        count() %>%
        mutate(labels=as.integer(is_driver) + 1,
               labels=c("non-driver", "driver")[labels]) %>%
        mutate(weight=ifelse(is_driver, w[2], w[1]),
               weight=weight/(sum(weight)),
               alpha=alpha,
               alpha2=alpha * weight/n,
               cutoff_nmodel=alpha2/n) %>%
        select(is_driver, cutoff_nmodel) %>%
        left_join(correct_ntests, by="is_driver")
    x3 <- x2 %>%
        mutate(is_driver=nonras %in% drivers) %>%
        left_join(cutoffs, by="is_driver") %>%
        mutate(cutoff_nmodel=round(-log10(cutoff_nmodel), 2),
               cutoff_ntest=round(-log10(cutoff_ntest), 2)) %>%
        mutate(p=round(-log10(p), 2),
               chisq_stat=round(chisq_stat, 3),
               chisq_p=round(-log10(chisq_p), 2),
               chisq_padj=round(-log10(chisq_padj), 2))
    x3
}

correct_overall <- function(overall, drivers, w, alpha){
    correct_ntests <- overall %>%
        mutate(is_driver=nonras %in% drivers) %>%
        group_by(is_driver) %>%
        count() %>%
        mutate(labels=as.integer(is_driver) + 1,
               labels=c("non-driver", "driver")[labels]) %>%
        mutate(weight=ifelse(is_driver, w[2], w[1]),
               weight=weight/(sum(weight)),
               alpha=alpha,
               alpha2=alpha * weight/n,
               cutoff_ntest=alpha2/n) %>%
        select(is_driver, cutoff_ntest)
    cutoffs <- overall %>%
        group_by(nonras) %>%
        summarize(is_driver=nonras %in% drivers,
                  .groups="drop") %>%
        distinct() %>%
        group_by(is_driver) %>%
        count() %>%
        mutate(labels=as.integer(is_driver) + 1,
               labels=c("non-driver", "driver")[labels]) %>%
        mutate(weight=ifelse(is_driver, w[2], w[1]),
               weight=weight/(sum(weight)),
               alpha=alpha,
               alpha2=alpha * weight/n,
               cutoff_nmodel=alpha2/n) %>%
        select(is_driver, cutoff_nmodel) %>%
        left_join(correct_ntests, by="is_driver")
    cutoffs
}

##
## For stratified analyses
##
p_adjust2 <- function(x, drivers, w=c(0.1, 0.9), alpha=0.05,
                      min_overall=5){ ## minimum marginal frequency
    MIN <- min_overall
    xx <- x %>%
        mutate(z=mean/sd,
               p=ifelse(z < 0, 2*pnorm(z), 2*(1-pnorm(z))),
               p=ifelse(p==0, min(p[p > 0]), p)) %>%
        filter(number_ras > MIN, number_nonras > MIN)
    if(nrow(xx) == 0) return(NULL)
    ct <- select(xx, c("00", "10", "01", "11", "n"))
    ct[is.na(ct)] <- 0
    ct$n <- ct$`00` + ct$`10` + ct$`01` + ct$`11`
    chisqfun <- function(x){
        xx <- matrix(x, ncol=2)
        stats <- tidy(chisq.test(xx))
    }
    tmp <- select(ct, -n)
    chisqstats <- apply(tmp, 1, chisqfun) %>%
        map_dfr(bind_rows) %>%
        select(statistic, p.value) %>%
        set_colnames(c("chisq_stat", "chisq_p"))
    ct2 <- bind_cols(ct, chisqstats)
    x2 <- select(xx, -c("00", "10", "01", "11", "n")) %>%
        bind_cols(ct2) %>%
        mutate(chisq_padj=p.adjust(chisq_p, method="BH"))
    ##
    ## MT adjustment based on number of genes not number of strata
    ##
    overall <- filter(x2, variable=="overall")
    cutoffs <- correct_overall(overall, drivers, w, alpha)
    x3 <- x2 %>%
        mutate(is_driver=nonras %in% drivers) %>%
        left_join(cutoffs, by="is_driver") %>%
        mutate(cutoff_nmodel=round(-log10(cutoff_nmodel), 2),
               cutoff_ntest=round(-log10(cutoff_ntest), 2)) %>%
        mutate(p=round(-log10(p), 2),
               chisq_stat=round(chisq_stat, 3),
               chisq_p=round(-log10(chisq_p), 2),
               chisq_padj=round(-log10(chisq_padj), 2))
    x3
}

harmonize_stats <- function(object, script, category, drivers){
    nonrasgenes <- nonras_genes(object)
    mt <- here("output", "co-mutation", script,
               paste0(nonrasgenes, ".rds")) %>%
        "["(file.exists(.)) %>%
        map_dfr(readRDS) %>%
        rename(group=category) %>%
        mutate(category=category) %>%
        group_by(group, ras) %>%
        nest()
    ##
    ## For the multicancer model, we don't really care about the
    ## high-level mu_beta[3] coefficient
    ##
    mt.overall <- filter(mt, group=="overall")
    mt.overall$multicancer <- mt.overall %>%
        pull(data) %>%
        map(p_adjust, drivers, exclude="mu_beta[3]")
    mt.overall <- mt.overall %>%
        select(-data)
    ##
    ## For sex, race, tumor type strata, we do care about the overall
    ##
    mt.confounders <- filter(mt, group != "overall")
    stats <- mt.confounders$data %>%
        map(p_adjust2, drivers)
    mt.confounders$strata <- stats
    mt.confounders  <- mt.confounders %>%
        select(-data)
    colnames(mt.overall) <- c("group", "ras", "stan")
    colnames(mt.confounders) <- c("group", "ras", "stan")
    mt <- bind_rows(mt.overall, mt.confounders)
    mt
}

harmonize_stats2 <- function(object, script, category, drivers){
    x <- object %>%
        separate(uid, c("ras", "nonras", "cancer"),
                 sep="~", remove=FALSE)
    nonrasgenes <- unique(x$nonras)
    mt <- here("output", "co-mutation", script,
               paste0(nonrasgenes, ".rds")) %>%
        "["(map_lgl(., file.exists)) %>%
        map_dfr(readRDS) %>%
        rename(group=category) %>%
        mutate(category=category) %>%
        group_by(group, ras) %>%
        nest()
    ##
    ## For TMB and mut sig., determine cutoff based on overall
    ##
    cat_  <- unique(mt$data[[1]]$category)
    if(cat_=="Rearrangement"){
        mt$stan <- mt$data %>%
            map(p_adjust, drivers)

    } else {
        mt$stan <- mt$data %>%
            map(p_adjust2, drivers)
    }
    mt <- select(mt, -data)
    mt
}

data_wrangle <- function(x){
    ## any mutation, tmb, and mutation signature should be a single category
    tmp <- x %>%
        filter(##cancer=="non-small cell lung cancer",
            group != "overall",
            variable != "overall")
    sigs <- c("1A", "1B", 2:21, "R2", "R3", "U1", "U2")
    vars <- tibble(variable=c("Female", "Male",
                              "Asign", "Black", "Unknown", "White",
                              "Primary", "Metastasis",
                              "2", "3", "4", "5", "6",
                              paste0("Signature.", sigs))) %>%
        mutate(variable2=case_when(variable=="2"~"TMB q2",
                                   variable=="3"~"TMB q3",
                                   variable=="4"~"TMB q4",
                                   variable=="5"~"TMB q5",
                                   variable=="6"~"TMB q6",
                                   variable=="Unknown"~"Race N/A",
                                   TRUE~variable),
               variable2=str_replace_all(variable2, "\\.", " "))
    rlevels <- tibble(ras=c("RAS_12-61", "RAS 12", "KRAS_12C"),
                      ras2=c("RAS Codons 12, 13, or 61",
                             "RAS Codon 12",
                             "KRAS G12C"))
    tmp <- filter(tmp, ras %in% rlevels$ras) %>%
        left_join(rlevels, by="ras") %>%
        mutate(ras=ras2,
               ras=factor(ras, rlevels$ras2)) %>%
        left_join(vars, by="variable") %>%
        filter(!is.na(variable2)) %>%
        mutate(variable2=factor(variable2, rev(vars$variable2))) %>%
        filter(category %in% c("Any mutation",
                               "Inactivating mutation"))
    tmp
}


gg_tmb <- function(x, xlimit){
    uid <- x[1, ] %>%
        unite(uid, c(ras, nonras, cancer), sep=",") %>%
        pull(uid)
    x %>%
        ggplot(aes(mean, label)) +
        geom_errorbarh(aes(xmin=`5%`, xmax=`95%`),
                       height=0.1, color="steelblue") +
        geom_point(shape=21, fill="white") +
        theme_bw(base_size=15) +
        xlim(xlimit) +
        ylab("") +
        xlab("") +
        ##xlab("Odds ratio") +
        ##facet_grid(cancer~`co-mutation`) +
        geom_vline(xintercept=1, linetype="dashed") +
        geom_text(aes(x=xlimit[2], y=label, label=n),
                  size=3, hjust=1) +
        theme(axis.text=element_text(size=8),
              strip.text=element_blank(),
              panel.grid=element_blank(),
              axis.text.y=element_text(lineheight=0.5),
              plot.title=element_text(size=6, lineheight=0.75)) +
        facet_grid(is_overall~., scales="free_y", space="free_y") +
        ggtitle(uid)
}
