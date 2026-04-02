
source("deepdive/00_a_load_data.R")
# yields:
# - count_data_synth
# - count_data_real
# - count_data_synth_long
# - count_data_real_long


library(tidyverse)
library(ggplot2)

################################################################################
# PCA on synthetic count data
# 
# PCA expects samples as rows, genes as columns, so we transpose.
# We also scale and center the data (z-scores), which is standard practice.

pca_result <- prcomp(t(count_data_synth), center = TRUE, scale. = TRUE)

# How much variance does each PC explain?
pca_var <- summary(pca_result)$importance
print(pca_var[, 1:4])

# Plot: variance explained per PC
df_var_explained <- data.frame(
    PC = paste0("PC", 1:ncol(pca_var)),
    variance_explained = pca_var["Proportion of Variance", ]
) |> mutate(PC = factor(PC, levels = PC))

p_scree <- ggplot(df_var_explained, aes(x = PC, y = variance_explained)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    ylab("Proportion of variance explained") +
    ggtitle("Scree plot - synthetic data")
p_scree

################################################################################
# Plot: PCA scatter (PC1 vs PC2)

df_scores <- as.data.frame(pca_result$x) |>
    mutate(sample = rownames(pca_result$x),
           condition = ifelse(str_detect(sample, "^WT"), "WT", "CompoundX"))

p_pca <- ggplot(df_scores, aes(x = PC1, y = PC2, color = condition)) +
    geom_point(size = 3) +
    theme_classic() +
    xlab(paste0("PC1 (", round(pca_var["Proportion of Variance", 1] * 100, 1), "%)")) +
    ylab(paste0("PC2 (", round(pca_var["Proportion of Variance", 2] * 100, 1), "%)")) +
    ggtitle("PCA - synthetic data")
p_pca

################################################################################
# Gene loadings: which genes contribute most to PC1 and PC2?
# The rotation matrix contains the loadings (weights) per gene per PC.

df_loadings <- as.data.frame(pca_result$rotation) |>
    rownames_to_column("gene") |>
    mutate(PC1_abs = abs(PC1), PC2_abs = abs(PC2))

# Top 15 genes by PC1 loading
top_pc1 <- df_loadings |> arrange(desc(PC1_abs)) |> head(15)
print("Top genes contributing to PC1:")
print(top_pc1 |> select(gene, PC1))

# Top 15 genes by PC2 loading
top_pc2 <- df_loadings |> arrange(desc(PC2_abs)) |> head(15)
print("Top genes contributing to PC2:")
print(top_pc2 |> select(gene, PC2))

# Plot loadings for PC1
p_loadings <- top_pc1 |>
    mutate(gene = fct_reorder(gene, PC1)) |>
    ggplot(aes(x = gene, y = PC1)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_classic() +
    ggtitle("Top 15 gene loadings on PC1")
p_loadings

################################################################################
# How is tpiA loaded across all PCs?

df_tpiA <- df_loadings |>
    filter(gene == "tpiA") |>
    select(gene, starts_with("PC"), -PC1_abs, -PC2_abs) |>
    pivot_longer(-gene, names_to = "PC", values_to = "loading") |>
    mutate(PC = factor(PC, levels = paste0("PC", 1:8)))

p_tpiA <- ggplot(df_tpiA, aes(x = PC, y = loading)) +
    geom_bar(stat = "identity") +
    theme_classic() +
    ggtitle("tpiA loading across all PCs")
p_tpiA

# Now also check out one of the noise genes

