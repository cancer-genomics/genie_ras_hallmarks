library(ggplot2)
library(ggbeeswarm)
#----------------------------------#
#' Adds count annotation to boxplot
#'
stat_box_data <- function(y, annot.level = -20) {
  return(
    data.frame(
      y = annot.level,
      label = paste0(length(y))
    )
  )
}

#-----------------------------------------#
# perform single feature coxph analysis,
# using top 1/3 of risk as one group, and bottom 2/3 as
# the other group
single.feature.km <- function(data, s, s.cens, param, title, cols = c("#377EB8", "#E41A1C")){
  data$surv = data[,s]
  data$surv.cens = data[,s.cens]
  data$var = data[,param]
  data <- subset(data, ! is.na(var))

  cox.fit <- coxph(Surv(surv,surv.cens) ~ var, data = data)
  riskscore <- cox.fit$linear.predictors
  riskcat <- ifelse(riskscore > quantile(riskscore, 0.667), 1, 0)

  if (length(unique(riskcat)) == 1){
    out <- list(plot = NULL, summary = paste(c(title, rep(NA, 4)), collapse = '\n '))
    return(out)
  }

  out.data <- data[,c(s, s.cens)]
  out.data[,param] <- riskcat
  plot <- km.plot(out.data, s, s.cens, param, title, cols)
  return(plot)
}
#-----------------------------------------#
# cap the survival field at user specified interval
# and adjust the censoring variable where necessary
cap.survival <- function(data, s, s.cens, threshold){
  fix.indices <- which(data[,s] > threshold)
  data[fix.indices, s] <- threshold
  data[fix.indices, s.cens] <- 0
  return(data)
}
#-----------------------------------------#
#' KM plots for survival analysis
#'
km.plot <- function(data, s, s.cens, param, title, cols){
  data$surv = data[,s]
  data$surv.cens = data[,s.cens]
  data$CAT = data[,param]

  cox.fit <- coxph(Surv(surv,surv.cens) ~ CAT, data = data)
  cox.pvalue = formatC(summary(cox.fit)$waldtest[3], digits = 4)

  diff.fit <- survdiff(Surv(surv, surv.cens) ~ CAT, data = data)
  logrank.pvalue = formatC(1 - pchisq(diff.fit$chisq, length(diff.fit$n) - 1), digits = 4)

  cox.survfit.km <- survfit(Surv(surv,surv.cens) ~ CAT, data = data)

  meds <- sapply(summary(cox.survfit.km)$table[,'median'],
                 function(x) sprintf('%3.2f',x))
  counts <- sapply(summary(cox.survfit.km)$table[,'records'],
                   function(x) sprintf('%d',x))


  title.text = paste(title, "\n",
                     "N =", paste(counts, collapse = ' vs '), "\n",
                     "P (Wald) =",cox.pvalue, "\n",
                     "P (Log-Rank) =", logrank.pvalue, "\n",
                     "Median Survivals = ", paste(meds, collapse = ' vs '))

  # print(sprintf('%s   -    %s', cox.pvalue, logrank.pvalue))

  if (length(unique(data$CAT)) == 2){
    title.text = paste0(title.text, '\n',
                        'Hazard Ratio: ', round(summary(cox.fit)$conf.int[1], digits = 2),
                        '(', round(summary(cox.fit)$conf.int[3], 2), ',',
                        round(summary(cox.fit)$conf.int[4], 2), ')')
  }

  p1 <- ggsurvplot(cox.survfit.km, data = data, risk.table = TRUE, palette = cols,
                   tables.height = 0.2, fontsize=12, tables.y.text = FALSE,
                   tables.theme = theme_cleantable(), risk.table.fontsize = 4.0,
                   risk.table.y.text.col = T, risk.table.y.text = FALSE,
                   font.legend = c(12, "plain", "black"),
                   ylab ="Survival (Probability)", xlab="Months",
                   title=title.text,
                   fontsize.title = 1, risk.table.col = "strata",
                   risk.table.title=NULL,
                   font.title = c(12, "plain", "black"),
                   conf.int=FALSE, break.time.by = 10,
                   font.x = c(12, "plain", "black"),
                   font.y = c(12, "plain", "black"),
                   font.xtickslab = c(12, "plain", "black"),
                   font.ytickslab = c(12, "plain", "black"))

  # return(p1)
  out <- list(plot = p1, summary = title.text,
              pval.wald= cox.pvalue,
              pval.logrank = logrank.pvalue)
  return(out)
}

#-----------------------------------------#
#' compare two variables in a cohort
#'
#' @import ggplot2
#' @import RColorBrewer
compare_feature_distribution <- function(data, variable, title = NULL ,  log.scale = FALSE, annot.level = 0.5){

  data <- data[,c('clinical.benefit', variable)]
  data <- na.omit(data)

  res <- data$clinical.benefit
  res <- factor(as.character(res), levels = c('DCB', 'NDB'))
  value <- data[,variable]

  f <- function(x) ifelse(x < 0.01, formatC(x, format = "e", digits = 2), round(x, 2))
  annot.lines <- c(paste0('Feature: ', variable),
  				   paste0('n = ', length(res)),
                   paste0('MW p-value = ', f(wilcox.test(value[res == 'DCB'], value[res == 'NDB'])$p.value)))


  if (! is.null(title)){
    annot.lines <- c(title, annot.lines)
  }


  plot.data <- data.frame(response = res,
                          value = value,
                          feature = variable)

  p <- ggplot(plot.data, aes(x = response, y = value, colour = response)) +
    geom_boxplot() +
    geom_quasirandom(width = 0.1) +
    theme_classic() +
    theme(legend.title = element_blank()) +
    scale_colour_brewer(palette = 'Set1') +
    labs(y = '', x= '',  title = paste(annot.lines, collapse = '\n'))

  if (log.scale == TRUE){
    p <- p + scale_y_log10()
  }

  p <- p + stat_summary(fun.data = function(x) stat_box_data(x, annot.level = annot.level),
                        geom = "text", hjust = 0.5, vjust = 1, color = 'black')
  return(p)
}
