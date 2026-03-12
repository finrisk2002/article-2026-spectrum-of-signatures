
library(mia)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(gridExtra)

nmf_models <- metadata(tseES)$NMF
variables <- as.data.frame(colData(tseES)) %>%
  select(c("DEATH", "DEATH_AGEDIFF", "DIAB", "DIAB_AGEDIFF", "LIVERDIS",
           "LIVERDIS_AGEDIFF", "AB1_SEPSIS_BACT", "AB1_SEPSIS_BACT_AGEDIFF"))

causes <- c("DEATH", "DIAB", "LIVERDIS", "AB1_SEPSIS_BACT")

all_results <- list()
assumptions <- c()
for (e in 1:length(causes)) {
  event <- causes[[e]]
  if (event == "DEATH") {
    model <- c(1:20)
  } else {
    model <- c(5,20)
  }
  results <- matrix(nrow = 0, ncol = 10)
  colnames(results) <- c("k", "name", "HR", "lower95", "upper95", "coeff", 
                         "p_val", "concordance", "ncomps", "q_val")
  
  for (i in 1:length(model)) {
    ncomps <- model[i]
    
    comps <- nmf_models$fit[[ncomps]]@fit@W
    taxa <- as.data.frame(nmf_models$fit[[ncomps]]@fit@H)
    
    # Name the components based on the maximum contribution
    rownames(taxa) <- colnames(taxa)[apply(taxa, 1, which.max)]
    colnames(comps) <- rownames(taxa)
    colnames(comps) <- gsub(" ", "_", colnames(comps))
    colnames(comps) <- gsub("-", "_", colnames(comps))
    
    all_covariates <- cbind(variables, comps)
    
    results_per_k <- matrix(nrow = 0, ncol = 9)
    colnames(results_per_k) <- c("k", "name", "HR", "lower95", "upper95", 
                                 "coeff", "p_val", "concordance", 
                                 "ncomps")
    
    all_covariates$event <- all_covariates[[causes[[e]]]]
    all_covariates$timediff <- all_covariates[[paste0(causes[[e]], "_AGEDIFF")]]
    all_covariates <- all_covariates %>% filter(timediff > 0)
    surv_response <- paste("Surv(timediff, event)")
    
    fml <- as.formula(paste(surv_response, "~", paste(colnames(comps), 
                                                      collapse = " + ")))
    smr <- summary(coxph(fml, data = all_covariates))
    # Get GLOBAL to evaluate if Cox prop haz assumptions are met
    assumptions <- c(assumptions, cox.zph(coxph(fml, 
                                   data = all_covariates))$table[i+1,3])
    for (k in 1:ncomps) {
      taxon <- rownames(smr$coefficients)[k]
      p_val <- smr$coefficients[[k,5]]
      coeff <- as.numeric(smr$coefficients[[k,1]])
      hr <- exp(coeff*1e6)
      se <- as.numeric(smr$coefficients[[k,3]])
      low95 <- exp(coeff*1e6 - 1.96*se*1e6)
      up95 <- exp(coeff*1e6 + 1.96*se*1e6)
      c_index <- smr$concordance

      add <- c(as.numeric(k), taxon, hr, low95, up95, coeff, as.numeric(p_val), 
               as.numeric(c_index[1]), ncomps)
      results_per_k <- rbind(results_per_k, add)
    }
    results_per_k <- as.data.frame(results_per_k)
    results_per_k$p_val <- as.numeric(results_per_k$p_val)
    results <- rbind(results, results_per_k)
  }
  results <- as.data.frame(results)
  results$coeff <- as.numeric(results$coeff)
  results$ncomps <- as.numeric(results$ncomps)
  results$concordance <- as.numeric(results$concordance)
  all_results[[event]] <- results
}
# 1-2 significant are expected with 26 tests (0.05*26=1.3)
sum(assumptions < 0.05)

# PROSPECTIVE RESULTS
df_end <- bind_rows(all_results, .id = "endpoint")
df_end <- df_end %>% mutate(endpoint = recode(endpoint,
                                              `DEATH` = "Death",
                                              `DIAB` = "Diabetes",
                                              `LIVERDIS` = "Liver dis.",
                                              `AB1_SEPSIS_BACT` = "Sepsis"))
df_end$q_val <-  p.adjust(df_end$p_val, method = "fdr")

write.csv(df_end, "./tables/survival.csv")

## Plotting
nk <- 20 # SELECT HERE EITHER 5 OR 20
df_end$plot <- es_names(df_end$name, nk)
df_end_nk <- df_end %>% filter(ncomps == nk)

if (nk == 5) {
  set.seed(15)
  b <- risk_plot(df_end_nk, "B", nk)
} else {
  set.seed(10)
  d <- risk_plot(df_end_nk, "D", nk)
}

# KM, categorical
comps <- nmf_models$fit[[nk]]@fit@W
taxa <- as.data.frame(nmf_models$fit[[nk]]@fit@H)
colnames(comps) <- colnames(taxa)[apply(taxa, 1, which.max)]
colnames(comps) <- es_names(colnames(comps), nk)

if (nk == 5) {
  variables$topic <- colnames(comps)[apply(comps, 1, which.max)]
  variables$topic <- as.factor(variables$topic)
  variables$topic <- relevel(variables$topic, ref ="ES-Bact")
  
  fit_nmf <- coxph(Surv(DEATH_AGEDIFF, DEATH) ~ topic, data = variables)
  s <- summary(fit_nmf)
  coef_table <- as.data.frame(s$coefficients)
  coef_table$q <- p.adjust(coef_table$`Pr(>|z|)`, "fdr")
  
  levels <- c("ES-Bact", rownames(coef_table)[coef_table$q < 0.05])
  levels <- gsub("topic", "", levels)
  
  newdata <- variables %>% filter(topic %in% c("ES-Bact", levels)) %>% 
    as.data.frame()
  newdata$topic <- case_when(
    newdata$topic == "ES-Esch" ~ paste0("ES-Esch, n=", 
                                          sum(newdata$topic=="ES-Esch")),
    newdata$topic == "ES-Prev" ~ paste0("ES-Prev, n=", 
                                          sum(newdata$topic=="ES-Prev")),
    newdata$topic == "ES-Bact" ~ paste0("Reference:\nES-Bact, n=", 
                                          sum(newdata$topic=="ES-Bact")))
  c <- km_plot(newdata, "C", nk)
  
  coef_table$low95 <- exp(coef_table$coef - 1.96*coef_table$`se(coef)`)
  coef_table$up95 <- exp(coef_table$coef + 1.96*coef_table$`se(coef)`)
  
} else {
  set.seed(3)
  variables$topic <- colnames(comps)[apply(comps, 1, which.max)]
  variables$topic <- as.factor(variables$topic)
  variables$topic <- relevel(variables$topic, ref ="ES-Bact. H uniformis")
  
  fit_nmf <- coxph(Surv(DEATH_AGEDIFF, DEATH) ~ topic, data = variables)
  s <- summary(fit_nmf)
  coef_table <- as.data.frame(s$coefficients)
  coef_table$q <- p.adjust(coef_table$`Pr(>|z|)`, "fdr")
  
  levels <- c("ES-Bact. H uniformis", rownames(coef_table)[coef_table$q < 0.05])
  levels <- gsub("topic", "", levels)
  
  newdata <- variables %>% 
    filter(topic %in% c("ES-Bact. H uniformis", levels)) %>% 
    as.data.frame()
  newdata$topic <- case_when(
    newdata$topic=="ES-Esch. ruysiae"~paste0("ES-Esch. ruysiae, n=",
                                           sum(newdata$topic=="ES-Esch. ruysiae")),
    newdata$topic=="ES-Prev. copri"~paste0("ES-Prev. copri, n=",
                                           sum(newdata$topic=="ES-Prev. copri")),
    newdata$topic=="ES-Prev. sp003447235"~paste0("ES-Prev. sp003447235, n=",
                                           sum(newdata$topic=="ES-Prev. sp003447235")),
    newdata$topic=="ES-Bact. H uniformis"~paste0("Reference:\nES-Bact. H uniformis, n=",
                                           sum(newdata$topic=="ES-Bact. H uniformis")),
    newdata$topic=="ES-Agat. faecis"~paste0("ES-Agat. faecis, n=",
                                                 sum(newdata$topic=="ES-Agat. faecis"))
    )
  e <- km_plot(newdata, "E", nk)
  
  coef_table$low95 <- exp(coef_table$coef - 1.96*coef_table$`se(coef)`)
  coef_table$up95 <- exp(coef_table$coef + 1.96*coef_table$`se(coef)`)
  
}
write.csv(coef_table,
          paste0("./tables/survival_categ_k", nk, ".csv"))

bd <- grid.arrange(b, d, ncol = 2, widths = c(0.45, 0.55))
ce <- grid.arrange(c, e, ncol = 2, widths = c(0.45, 0.55))
bcde <- grid.arrange(bd, ce, ncol = 1, heights = c(0.48, 0.52))

fig2 <- grid.arrange(fig2a, bcde, ncol = 1, heights = c(0.32, 0.68))
ggsave(paste0(output_dir, "fig2.png"), fig2, width = 9.5, height = 7.2,
       dpi = 300, bg = "white")
