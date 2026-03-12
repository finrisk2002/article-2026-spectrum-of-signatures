
library(mia)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# This script includes
# 1) Selection of k using cosine similarity
# 2) Overview plot (Figure 1)

tseES <- readRDS("TSE_GG2_MGS.rds")

# Selection of k using cosine similarity index
nmf_models <- metadata(tseES)$NMF
data <- assay(altExp(tseES, "prevalent"), "counts")

cosine_colwise <- function(X, Y) {
  numerators <- sum(X * Y)
  x_norms <- sqrt(sum(X*X))
  y_norms <- sqrt(sum(Y*Y))
  denominators <- x_norms * y_norms
  numerators / denominators}
cos_per_k <- c()
for (k in 1:20) {
  comps <- nmf_models$fit[[k]]@fit@W
  taxa <- nmf_models$fit[[k]]@fit@H
  reconstructed <- t(comps %*% taxa) # get reconstructed data matrix: counts
  cos_per_k[k] <- cosine_colwise(data, reconstructed)
}
df <- data.frame(index = seq_along(cos_per_k), value = cos_per_k)
df$index[df$value==max(df$value)]

# Overview plot (Figure 1)
nk <- 3
nmf_models <- metadata(tseES)$NMF
comps <- nmf_models$fit[[nk]]@fit@W
comps <- comps/rowSums(comps)

set.seed(50)
random_indices <- sample(1:nrow(comps), 8, replace = FALSE)

# Extract the sampled rows
rowsample <- comps[random_indices, ]
rowsample <- rowsample/rowSums(rowsample)
taxa <- as.data.frame(nmf_models$fit[[nk]]@fit@H)
rownames(taxa) <- c("ES-1", "ES-2", "ES-3")
colnames(rowsample) <- rownames(taxa)

breaks <- seq(0, 1, length.out = 5) # to make the colors more interpretable
breaks <- unique(c(0, breaks, 1))^3
col_fun <- colorRamp2(
  breaks,
  colorRampPalette(c("honeydew", "darkseagreen1", "darkseagreen2" , 
                     "darkseagreen3", "palegreen4"))(length(breaks))
)
h1 <- Heatmap(rowsample,
        name = "Memb.",
        show_heatmap_legend = FALSE,
        col = col_fun, 
        cluster_rows = FALSE, 
        show_row_dend = FALSE,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE, 
        column_names_rot = 0,
        column_names_gp = gpar(fontsize = 11),
        row_names_gp = gpar(fontsize = 11),
        column_title = "") # Sample membership
features <- c()
for (i in 1:nk) {
  features <- c(features, get_key(as.data.frame(t(taxa)), i, 2))
}
features <- unique(features)
coltaxa <- taxa[,features]
coltaxa <- as.matrix(coltaxa)
coltaxa <- coltaxa/rowSums(coltaxa)
h2 <- Heatmap(coltaxa, name = "Weight",
              col = col_fun, 
              cluster_rows = FALSE, 
              show_row_dend = FALSE,
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              show_row_names = FALSE,
              show_column_names = FALSE, 
              column_order = sort(colnames(coltaxa)),
              row_names_gp = gpar(fontsize = 11),
              heatmap_legend_param = list(
                at = c(min(coltaxa), max(coltaxa)),
                labels = c("Low", "High")),
              column_title = "")
abundances <- assay(altExp(tseES, "prevalent"), "relabundance")[, random_indices]
abundances <- t(subset(abundances, rownames(abundances) %in% features))
abundances <- abundances/rowSums(abundances)
breaks <- seq(0, 1, length.out = 5)
breaks <- unique(c(0, breaks, 1))^5
h0 <- Heatmap(abundances,
              name = "Abund.",
              show_heatmap_legend = TRUE,
              col = col_fun, 
              cluster_rows = FALSE, 
              show_row_dend = FALSE,
              cluster_columns = FALSE,
              show_column_dend = FALSE,
              column_order = sort(colnames(coltaxa)),
              show_row_names = FALSE,
              show_column_names = FALSE, 
              column_names_rot = 45,
              column_names_gp = gpar(fontsize = 11),
              row_names_gp = gpar(fontsize = 11),
              heatmap_legend_param = list(
                at = c(min(abundances), max(abundances)),
                labels = c("Low", "High")),
              column_title = "")
g0 <- grid.grabExpr(draw(h0, padding = unit(c(5, 10, 5, 30), "mm"), newpage = FALSE))
g1 <- grid.grabExpr(draw(h1, padding = unit(c(5, 10, 5, 15), "mm"), newpage = FALSE))
g2 <- grid.grabExpr(draw(h2, padding = unit(c(5, 10, 5, 5), "mm"), newpage = FALSE))

blank <- grid::nullGrob()  # empty grob
g3 <- grid.arrange(blank, g2, blank, heights = c(0.2, 0.45, 0.35), ncol = 1)

fig1 <- grid.arrange(g0, g1, g3, widths = c(0.45, 0.22, 0.36), ncol = 3)
ggsave(paste0(output_dir, "fig1_overview.png"), fig1, width = 9.5, height = 2.8,
       dpi = 300, bg = "white")
