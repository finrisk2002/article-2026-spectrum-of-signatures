
get_key <- function(df, index, n) {
  df$species <- rownames(df)
  drivers <- df[order(df[[index]], decreasing = TRUE),]$species[1:n]
  return(drivers)
}

es_names <- function(names, k) {
  if (k == 5) {
    names <- gsub(" ", "_", names)
    es <- case_when(
      names == "Agathobacter_rectalis" ~ "ES-Agat",
      names == "Bacteroides_H_cellulosilyticus"~ "ES-Esch",
      names == "Gemmiger_A_73129_qucibialis" ~ "ES-Alis",
      names == "Prevotella_copri" ~ "ES-Prev",
      names == "Bacteroides_H_uniformis" ~ "ES-Bact"
    )
  } else { # k == 20
    es <- gsub("_", " ", names)
    es <- gsub("\\b\\d+\\b", "", es)
    es <- sub("^(.{1,4})(?:[^ ]*)(\\s.*)$", "\\1.\\2", es,
              perl = TRUE)
    es <- gsub("  ", " ", es)
    es <- paste0("ES-", es)
  }
  return(es)
}

## FIG 2b,d
risk_plot <- function(df_end, title, k) {
  if (k == 5) {
    cols <- c("ES-Alis" = "steelblue2", "ES-Agat" = "turquoise3", "ES-Prev" = "indianred1", 
              "ES-Esch" = "gray20")
  } else {
    cols <- c("ES-Alis. A putredinis" = "steelblue2", "ES-Prev. copri" = "indianred1", 
              "ES-Prev. sp003447235" = "orange",
              "ES-Esch. ruysiae" = "gray20", "ES-Flav. plautii" = "deeppink", 
              "ES-RUG1. sp900066395" = "seagreen", 
              "ES-Phoc. A vulgatus" = "mediumpurple3")
  }
  nonsig <- df_end %>% filter(q_val > 0.05)
  nonsig$plot <- "qval>0.05" 
  subs_higher <- df_end %>% filter(coeff > 0 & q_val < 0.05)
  subs_lower <- df_end %>% filter(coeff < 0 & q_val < 0.05)
  
  breaks <- as.data.frame(seq(1.5, 3.5, by = 1))
  colnames(breaks) <- c("breaks")
  
  p <- ggplot() +
    geom_vline(data = breaks, 
               aes(xintercept = breaks), colour = "grey80") +
    geom_hline(yintercept=0.05, color = "grey40", linewidth=1) +
    # geom_jitter(nonsig, mapping = aes(x = endpoint, y = q_val, fill = plot),
    #             shape = 23, alpha = 0.7, size = 2.5, width = 0.3, height = 0) +
    geom_jitter(subs_lower, mapping = aes(x = endpoint, y = q_val, fill = plot),
                shape = 25, alpha = 0.8, size = 3.5, width = 0.3, height = 0) +
    geom_jitter(subs_higher, mapping = aes(x = endpoint, y = q_val, fill = plot),
                shape = 24, alpha = 0.8, size = 3.5, width = 0.3, height = 0) +
    scale_fill_manual(values = cols) +
    scale_y_log10(limits = c(1e-8, 1)) + 
    coord_cartesian(clip = "off") +
    theme_light(base_size = 14) +
    theme(plot.title = element_text(face = "bold", size = 24), axis.text.x = 
            element_text(angle = 20, hjust = 1), plot.margin = margin(10, 20, 10, 10),
          axis.title.y = element_text(vjust=3), legend.title = element_text(size = 12),
          legend.text = element_text(face = "italic")) +
    labs(x = "",
         y = "q-val",
         fill = "\u25B2HR>1;\u25BCHR<1",
         #size = "concordance",
         title = title,
         subtitle = paste0("k = ", k))
  return(p)
}
# Fig 2 c, e
km_plot <- function(newdata, title, k) {
  if (k == 5) {
    newdata$topic <- factor(newdata$topic, levels = c("Reference:\nES-Bact, n=1794", 
                                                      "ES-Prev, n=832",
                                                      "ES-Esch, n=807"))
    cols <- c("ES-Prev, n=832" = "indianred1", 
              "ES-Esch, n=807" = "gray20", "Reference:\nES-Bact, n=1794" = "violet")
  } else {
    newdata$topic <- factor(newdata$topic, levels = c("Reference:\nES-Bact. H uniformis, n=631",
                                                      "ES-Agat. faecis, n=465",
                                                      "ES-Prev. copri, n=578",
                                                      "ES-Prev. sp003447235, n=335",
                                                      "ES-Esch. ruysiae, n=161"))
    cols <- c("ES-Prev. copri, n=578" = "indianred1", 
              "ES-Prev. sp003447235, n=335" = "orange",
              "ES-Esch. ruysiae, n=161" = "gray20",
              "ES-Agat. faecis, n=465" = "chartreuse3",
              "Reference:\nES-Bact. H uniformis, n=631" = "violet")
  }
  p <- survfit2(Surv(DEATH_AGEDIFF, DEATH) ~ topic, data = newdata) |>
    ggsurvfit() +
    theme_light(base_size = 14) +
    add_confidence_interval() +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ylim(0.59, 1) +
    coord_cartesian(clip = "off") +
    theme(plot.title = element_text(face = "bold", size = 24),
          plot.margin = margin(10, 20, 10, 10), 
          legend.text = element_text(face = "italic")) +
    labs(x = "Years", y = "Survival probability", title = title,
         subtitle = paste0("k = ", k))
  return(p)
}

# inverse rank normalization
inv <- function(x) {
  r <- rank(x, na.last = "keep", ties.method = "average")
  x_rinv <- qnorm((r - 0.5) / sum(!is.na(x)))
  scale(x_rinv)
}