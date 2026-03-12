
library(mia)
library(dplyr)
library(ggplot2)
library(alto)
library(purrr)
library(stringr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(gridExtra)

nmf_models <- metadata(tseES)$NMF
all_models <- c(nmf_models$fit)

models <- lapply(seq_along(all_models), function(i) {
  res <- all_models[[i]]
  list(gamma = res@fit@W/rowSums(res@fit@W), 
       beta = 10*res@fit@H) # H to the same scale as topics from alto (10x)
})
names(models) <- c(seq_along(nmf_models$fit))

# USING ALTO
# Align the topics (components) across the NMF models
alignment <- align_topics(models) # method = product

source("../alto_script.R") # edited alto package

f3a <- plot_beta_edit(alignment, models = c(5, 20), filter_by = "beta", 
                      n_features = 60) + labs(title = "A") +
  theme(plot.title= element_text(face = "bold", size = 18), 
        axis.text.x = element_blank())

colnames(alignment@models$`20`$beta) <- es_names(colnames(alignment@models$`20`$beta), 20)
colnames(alignment@models$`20`$beta)[duplicated(colnames(alignment@models$`20`$beta))] <- c(1:8)

# Get prospective results from prospective.R
results <- df_end %>% filter(df_end$endpoint == "Death") 
signif_results <- results %>% filter(q_val < 0.05)

# locally edited functions from alto package
f3b <- plot_alignment_edit(alignment, add_leaves = TRUE, n_features_in_leaves = 1, 
                           leaves_text_size = 9, add_sig = signif_results) + 
  theme_light() + theme(text = element_text(size=9), 
                        plot.title = element_text(face = "bold", size = 18)) + 
  ggtitle("B") + # All-cause mortality
  annotate("rect", xmin = 19.5, xmax = 20.5, ymin = 0, ymax = 2,
           alpha = 0, color = "gray30", lwd = .5) + 
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = 0, ymax = 2,
           alpha = 0, color = "gray30", lwd = .5) +
  annotate("label", x = 6, y = 0, size = 9/.pt,
           label = "*q-val<0.05, **q-val<0.01, ***q-val<0.001")

# Similarity between k=5 and k=20
m1 <- nmf_models$fit[[5]]@fit@W
m2 <- nmf_models$fit[[20]]@fit@W

taxa <- as.data.frame(nmf_models$fit[[5]]@fit@H)
colnames(m1) <- colnames(taxa)[apply(taxa, 1, which.max)]
colnames(m1) <- es_names(colnames(m1), 5)

taxa2 <- as.data.frame(nmf_models$fit[[20]]@fit@H)
colnames(m2) <- colnames(taxa2)[apply(taxa2, 1, which.max)]
colnames(m2) <- es_names(colnames(m2), 20)

products <- t(m1) %*% m2
products <- products/max(products)

# plot new heatmap
breaks <- seq(0, 1, length.out = 5)
breaks <- unique(c(0, breaks, 1))
col_fun <- colorRamp2(
  breaks,
  colorRampPalette(c("snow", "thistle2", "thistle3", "thistle4"))(length(breaks))
)
f3c <- Heatmap(t(products), 
               name = "Rel. similarity",
               col = col_fun,
               cluster_rows = FALSE,
               show_row_dend = FALSE,
               cluster_columns = FALSE,
               show_column_dend = FALSE,
               show_row_names = TRUE,
               show_column_names = TRUE,
               column_names_rot = 0,
               column_names_gp = gpar(fontsize = 9),
               row_names_gp = gpar(fontsize = 9),
               column_title = "k=5",
               row_title = "k=20",
               row_title_gp = grid::gpar(fontface = "bold", fontsize = 11),
               column_title_gp = grid::gpar(fontface = "bold", fontsize = 11),
               column_names_centered = TRUE)

c <- grid.grabExpr(draw(f3c, padding = unit(c(5, 10, 5, 10), "mm"), 
                        column_title = "C",
                        column_title_gp = gpar(fontsize = 18, 
                                               fontface = "bold")))

right <- grid.arrange(f3b, c, ncol = 1, heights = c(0.52, 0.48))

fig3 <- grid.arrange(f3a, right, ncol = 2, widths = c(0.53, 0.47))
ggsave(paste0(output_dir, "fig3.png"), fig3, width = 14, height = 6.8,
       dpi = 300, bg = "white")
