
library(mia)
library(dplyr)
library(ggplot2)
library(rstatix)

### Visualizing FINRISK data in terms of ES
# age, sex, BMI, antibiotics, pregnancy, vegetarian diet
# prior event: diabetes, sepsis, liver diseases

# Kendall's tau index
nk <- 5
nmf_models <- metadata(tseES)$NMF
comps <- nmf_models$fit[[nk]]@fit@W
taxa <- as.data.frame(nmf_models$fit[[nk]]@fit@H)
rownames(taxa) <- colnames(taxa)[apply(taxa, 1, which.max)]
rownames(taxa) <- es_names(rownames(taxa), 5)
colnames(comps) <- rownames(taxa)

variables <- as.data.frame(colData(altExp(tseES, "prevalent"))) %>%
  select(c("BL_AGE", "MEN", "BMI", "FR02_100G", "BL_USE_RX_J01", "GRAVID",
           "DEATH", "DEATH_AGEDIFF","PREVAL_DIAB", "PREVAL_LIVERDIS", 
           "PREVAL_AB1_SEPSIS_BACT"))
variables$VEGE <- case_when(
  variables$FR02_100G == 1 ~ 0,
  variables$FR02_100G == 2 ~ 1,
  variables$FR02_100G == NA ~ NA)

variables$PREG <- case_when(
  variables$GRAVID == 1 ~ 0,
  variables$GRAVID == 2 ~ 1,
  variables$GRAVID == NA ~ NA)
  
data <- cbind(comps, variables)

vars <- c("BL_AGE", "MEN", "BMI", "VEGE", "PREVAL_DIAB", "PREVAL_LIVERDIS",
          "PREVAL_AB1_SEPSIS_BACT", "PREG", "BL_USE_RX_J01")
kendall <- matrix(nrow = 0, ncol = 5)

kendall <- data %>%
  cor_test(vars = starts_with("ES"), vars2 = all_of(vars),
           method = "kendall", pairwise.complete.obs = TRUE)

kendall$var2 <- case_when(
  kendall$var2 == "BL_AGE" ~ "Age",
  kendall$var2 == "MEN" ~ "Man",
  kendall$var2 == "BMI"~ "BMI",
  kendall$var2 == "BL_USE_RX_J01" ~ "Antibiotics (6 months)",
  kendall$var2 == "PREG" ~ "Pregnant",
  kendall$var2 == "VEGE" ~ "Vegetarian",
  kendall$var2 == "PREVAL_DIAB"~ "Prevalent diabetes",
  kendall$var2 == "PREVAL_LIVERDIS" ~ "Prior liver disease",
  kendall$var2 == "PREVAL_AB1_SEPSIS_BACT" ~ "Prior sepsis")

kendall <- kendall %>%
  mutate(sig_label = case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE ~ ""))

write.csv(kendall, "./tables/kendall.csv")

kendall <- kendall %>% 
  mutate(sig = case_when(
    p < 0.05 ~ "p-val<0.05",
    p >= 0.05 ~ "Not significant"))
kendall <- kendall %>% filter(!var2 %in% c("Vegetarian", "Pregnant")) # for all es non-significant

# 1100, 250
fig2a <- ggplot(kendall, aes(x = var2, y = cor, fill = sig)) +
  geom_bar(stat = "identity") +
  ylim(c(-0.15, 0.15)) +
  scale_fill_manual(
    values = c(
      "p-val<0.05"     = "hotpink3",
      "Not significant" = "azure4"
    ),
    name = "") +
  facet_wrap(~ var1, scales = "free_x", nrow = 1) +
  coord_flip() +
  theme_bw(base_size = 14) +
  theme(legend.position = "right",
        plot.title = element_text(face = "bold", size = 24)) +
  labs(x = "", y = "Kendall's tau", title = "A", fill = "")
