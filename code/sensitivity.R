
library(mia)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(gridExtra)

nmf_models <- metadata(tseES)$NMF
variables <- as.data.frame(colData(tseES)) %>%
  select(c("DEATH", "DEATH_AGEDIFF"))
altExp(tseES, "genus") <- agglomerateByPrevalence(tseES, rank="genus", 
                                                  assay.type="relabundance",
                                                  detection = 0.01/100, 
                                                  prevalence = 10/100)
nmf20 <- nmf_models$fit[[20]]@fit@W
taxa <- nmf_models$fit[[20]]@fit@H
colnames(nmf20)  <- colnames(taxa)[apply(taxa, 1, which.max)]

names <- unique(colnames(nmf20))
taxainfo <- as.data.frame(rowData(altExp(tseES, "prevalent"))@listData)
gen <- unique(taxainfo$genus[taxainfo$species %in% names])

df_genus <- assay(altExp(tseES, "genus"), "relabundance")
df_genus <- t(df_genus[gen,])
df_es <- nmf20

all_results <- list()
assumptions <- c()
for (e in 1:2) {
  event <-"DEATH"
  
  if (e == 1) {
    df_sens <- df_es
  } else {
    df_sens <- df_genus
  }
  colnames(df_sens) <- gsub("-", "?", colnames(df_sens))
  colnames(df_sens) <- gsub(" ", "_", colnames(df_sens))
  
  df_norm <- df_sens
  # # Inverse rank transformation
  for (j in 1:(ncol(df_sens))) df_norm[, j] <- inv(df_sens[, j])
  colVars(df_norm)
  all_covariates <- cbind(variables, df_norm)
  
  surv_response <- paste("Surv(DEATH_AGEDIFF, DEATH)")
  results <- matrix(nrow = 0, ncol = 7)
  
  for (i in 1:ncol(df_sens)) {
    predictor <- colnames(df_sens)[i]
    fml <- as.formula(paste(surv_response, "~", predictor))
    smr <- summary(coxph(fml, data = all_covariates))
    # Get GLOBAL to evaluate if Cox prop haz assumptions are met
    # assumptions <- c(assumptions, cox.zph(coxph(fml,
    #                                data = all_covariates))$table[1,3])
    # if (cox.zph(coxph(fml, data = all_covariates))$table[1,3]<0.05) {
    #   print(plot(cox.zph(coxph(fml, data = all_covariates))))
    # }
    p_val <- as.numeric(smr$coefficients[[1,5]])
    hr <- as.numeric(smr$coefficients[[1,2]])
    se <- as.numeric(smr$coefficients[[1,3]])
    low95 <- exp(as.numeric(smr$coefficients[[1,1]])-1.96*as.numeric(smr$coefficients[[1,3]]))
    up95 <- exp(as.numeric(smr$coefficients[[1,1]])+1.96*as.numeric(smr$coefficients[[1,3]]))
    concordance <- smr$concordance[1]

    add <- c(predictor, hr, low95, up95, se, p_val, concordance)
    results <- rbind(results, add)
  }
  results <- as.data.frame(results)
  colnames(results) <- c("name", "HR", "lower95", "upper95", "se", 
                         "p_val", "concordance")
  results$q_val <- p.adjust(results$p_val, "fdr")
  if (e == 1) {
    results$name <- es_names(results$name, 20)
    all_results[["ES"]] <- results
  } else {
    results$name <- gsub("_", " ", results$name)
    all_results[["Genus"]] <- results
  }
}
sum(assumptions < 0.05) # OK


# PROSPECTIVE RESULTS
df <- bind_rows(all_results, .id = "Set")
df$HR <- as.numeric(df$HR)
df$significant <- case_when(df$q_val < 0.05 & df$q_val > 0.01 ~ "*",
                            df$q_val < 0.01 & df$q_val > 0.001 ~ "**",
                            df$q_val < 0.001 ~ "***",
                            df$q_val > 0.05 ~ "")

write.csv(df, "./tables/univariate_sensitivity.csv")
