
library(mia)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(gridExtra)

# This script includes
# 1) Comparison of the effect sizes to mortality between ESs and their representative genera
# 2) Adjusting the model effects for the covariates in FR02 subset

# Sensitivity analysis 1: compare the effect sizes of genera and ESs
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
write.csv(df, "./tables/univariate.csv")


# Sensitivity analysis 2: adjust for covariates,subset FR02
# k=5 and k=20
nmf_models <- metadata(tseES)$NMF
variables <- as.data.frame(colData(tseES)) %>%
  select(c("DEATH", "DEATH_AGEDIFF", "BL_AGE", "MEN", "BMI", "BL_USE_RX_J01",
           "FID", "Barcode"))
variables$FR02 <- as.factor(ifelse(startsWith(variables$FID, "40"),
                             0, 1))
variables <- subset(variables, variables$FR02 == 1)

variables$MEN <- as.factor(variables$MEN)

all_results <- list()
assumptions <- c()
adjust <- c("FR02", "Covariates")

for (e in adjust) {
  event <- "DEATH"
  model <- c(5,20)

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
    
    comps <- subset(comps, rownames(comps) %in% variables$Barcode)
    all_covariates <- cbind(variables, comps)
    
    results_per_k <- matrix(nrow = 0, ncol = 9)
    colnames(results_per_k) <- c("k", "name", "HR", "lower95", "upper95", 
                                 "coeff", "p_val", "concordance", 
                                 "ncomps")
    
    all_covariates$event <- variables$DEATH
    all_covariates$timediff <- variables$DEATH_AGEDIFF
    surv_response <- paste("Surv(timediff, event)")
    
    if (e == "Covariates") {
      fml <- as.formula(paste(surv_response, "~", paste(colnames(comps), 
                                                        collapse = " + "),
                              "+BL_AGE+MEN+BMI+BL_USE_RX_J01"))
      ncovs <- 4
    } else {
      fml <- as.formula(paste(surv_response, "~", paste(colnames(comps), 
                                                        collapse = " + ")))
      ncovs <- 0
    }
    smr <- summary(coxph(fml, data = all_covariates))
    for (k in 1:(ncomps+ncovs)) {
      pred <- rownames(smr$coefficients)[k]
      p_val <- smr$coefficients[[k,5]]
      coeff <- as.numeric(smr$coefficients[[k,1]])
      
      if (k < ncomps) { # ESs
      hr <- exp(coeff*1e6)
      se <- as.numeric(smr$coefficients[[k,3]])
      low95 <- exp(coeff*1e6 - 1.96*se*1e6)
      up95 <- exp(coeff*1e6 + 1.96*se*1e6)
      } else { # covariates
        hr <- exp(coeff)
        se <- as.numeric(smr$coefficients[[k,3]])
        low95 <- exp(coeff - 1.96*se)
        up95 <- exp(coeff + 1.96*se)
      }
      c_index <- smr$concordance
      
      add <- c(as.numeric(k), pred, hr, low95, up95, coeff, as.numeric(p_val), 
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
  all_results[[e]] <- results
}

# RESULTS
df_end <- bind_rows(all_results, .id = "endpoint")
df_end <- df_end %>% mutate(endpoint = recode(endpoint,
                                              `FR02` = "Study year",
                                              `Covariates` = "Covariates"))
df_end$HR <- as.numeric(df_end$HR)
df_end$lower95 <- as.numeric(df_end$lower95)
df_end$upper95 <- as.numeric(df_end$upper95)

df_end$name <- case_when(
  df_end$name == "BL_AGE" ~ "Age",
  df_end$name == "MEN1" ~ "Man",
  df_end$name == "BMI"~ "BMI",
  df_end$name == "BL_USE_RX_J01" ~ "Antibiotics (6 months)",
  TRUE ~ df_end$name)

df_end$name[1:5] <- es_names(df_end$name[1:5], 5)
df_end$name[6:25] <- es_names(df_end$name[6:25], 20)
df_end$name[26:30] <- es_names(df_end$name[26:30], 5)
df_end$name[35:54] <- es_names(df_end$name[35:54], 20)

# supplement figure
subs1 <- subset(df_end, ncomps == 5 & endpoint == "Study year")
s1 <- ggplot(subs1, aes(x = HR, y = name)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_light(base_size = 14) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", y = "",
    title = "A", subtitle = "FR2002 cohort, k=5") + 
  theme(axis.text.y = element_text(face = "italic"),
                         plot.title = element_text(face = "bold", size = 24),
        plot.margin = unit(c(10, 10, 30, 90), "pt"))
subs2 <- subset(df_end, ncomps == 20 & endpoint == "Study year")
s2 <- ggplot(subs2, aes(x = HR, y = name)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_light(base_size = 14) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", y = "",
       title = "B", subtitle = "FR2002 cohort, k=20") + 
  theme(axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 24),
        plot.margin = unit(c(10, 10, 60, 10), "pt"))
subs3 <- subset(df_end, ncomps == 5 & endpoint == "Covariates")
s3 <- ggplot(subs3, aes(x = HR, y = name)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_light(base_size = 14) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", y = "",
       title = "C", subtitle = "+ covariates") + 
  theme(axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 24),
        plot.margin = unit(c(10, 10, 10, 30), "pt"))
subs4 <- subset(df_end, ncomps == 20 & endpoint == "Covariates")
s4 <- ggplot(subs4, aes(x = HR, y = name)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_light(base_size = 14) +
  scale_x_log10() +
  labs(x = "Hazard Ratio (log scale)", y = "",
       title = "D", subtitle = "+ covariates") + 
  theme(axis.text.y = element_text(face = "italic"),
        plot.title = element_text(face = "bold", size = 24))

left <- grid.arrange(s1, s2, heights = c(0.35, 0.65), ncol = 1)
right <- grid.arrange(s3, s4, heights = c(0.35, 0.65), ncol = 1)
sens <- grid.arrange(left, right, ncol = 2, widths = c(0.45, 0.55))

ggsave(paste0(output_dir, "sensitivity_fr02.png"), sens, width = 14, height = 10,
       dpi = 300, bg = "white")
