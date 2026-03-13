
library(NMF)
library(mia)

tse <- readRDS("TSE_GG2_MGS.RDS")

# filter out the samples with read count less than 50 000
df <- as.data.frame(assay(tse, "counts"))
total_reads <- colSums(assay(tse, "counts"))
to_prune <- total_reads > 50000
filtered <- colData(tse)[to_prune,]
tse <- tse[,rownames(filtered)]

# filter out the least prevalent features
altExp(tse, "prevalent") <- agglomerateByPrevalence(tse, rank = "species",
                                                    assay.type = "relabundance",
                                                    prevalence = 10/100,
                                                    detection = 0.01/100)
# Filter out the "Other" taxa
altExp(tse, "prevalent") <- altExp(tse, 
                                   "prevalent")[!(rownames(altExp(tse, 
                                               "prevalent")) %in% c("Other")),]
# NMF
vec <- c(1:20) # k
x <- t(assay(altExp(tse, "prevalent"), "counts")) # species
# nmf across k
models <- nmf(x, rank=vec, nrun = 8, seed = 12345)

metadata(tse)$NMF <- models

