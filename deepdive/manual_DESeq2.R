################################################################################

source("deepdive/import_data.R")
# - count_data_synth
# - count_data_real
# - count_data_synth_long
# - count_data_real_long
# - df_dispersion_synth
# - df_dispersion_real

library(tidyverse)
library(ggplot2)

library(DESeq2)

################################################################################
# Let's also do some manual DESeq2-similar things
# First determine total reads per sample
df_total_reads <- count_data_real_long |>
    group_by(sample) |>
    summarize(total_counts = sum(counts))
    
# Now plot in bar graph
p <- df_total_reads |>
  mutate(sample=str_remove(sample, "Atha_WT_Col_0_sl_")) |>
  ggplot(aes(x=sample, y=total_counts)) +
    geom_bar(stat="identity") + 
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
# p    
# save
ggsave("plots/total_counts_per_sample.pdf", plot=p, width=7, height=7, units='cm')


################################################################################
# Now normalize them in DESeq2's way

# Estimate the "size factors"
# Create reference sample; geometric mean
df_size_factors <- count_data_real_long |>  
  # Calculate geometric means of gene values
  group_by(gene) |>
  mutate(geom_mean = median(prod(counts)^(1/(length(counts))))) |>
  ungroup() |>
  # Now normalize the counts by geom mean
  mutate(counts_norm = counts/geom_mean) |>
  # And calculate median normalized counts (=size factor)
  group_by(sample) |>
  summarize(size_factor = median(counts_norm, na.rm=TRUE))

# now plot those
p <- df_size_factors |>
  mutate(sample=str_remove(sample, "Atha_WT_Col_0_sl_")) |>
  ggplot(aes(x=sample, y=size_factor)) +
  geom_bar(stat="identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# p    
# save
ggsave("plots/size_factors_per_sample.pdf", plot=p, width=7, height=7, units='cm')


################################################################################
# Let's see how well this matches DESeq2's assessment
dds <- DESeqDataSetFromMatrix(countData = count_data_real,
                              colData   = count_data_real_meta,
                              design    = ~ condition)
dds <- estimateSizeFactors(dds)
the_sizefactors <- sizeFactors(dds)
df_size_factors_deseq2 <- data.frame(
  size_factor = unname(the_sizefactors),
  sample = names(the_sizefactors))

# now comparison
df_size_factors_all <-  
  rbind( df_size_factors |> mutate(origin="manual"), 
        df_size_factors_deseq2 |> mutate(origin = "deseq2"))
# now plot
p <- df_size_factors_all |>
  mutate(sample = str_remove(sample, "Atha_WT_Col_0_sl_")) |>
  ggplot(aes(x = sample, y = size_factor, fill = origin)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = c("manual" = "black", "deseq2" = "red"))
# p
ggsave("plots/size_factors_manual_vs_deseq2.pdf", plot = p, width = 9, height = 7, units = "cm")

# We can now also get normalized count matrix from deseq2
normalized_counts <- counts(dds, normalized=TRUE)
# Which doesn't have all sample totals equal
colSums(normalized_counts)

################################################################################
# PART II, PCA normalization
# See `deepdive/pca_example.R` for pca fundamentals

# There tends to be a trend in the variance, compared to mean expression values.
p <- ggplot(df_dispersion_real, aes(x=mean, y=variance)) +
    geom_point(color='grey') + 
    scale_x_log10() + scale_y_log10() + 
    theme_classic() +
    xlab('Mean gene expression (all conditions)') +
    ylab('Variance in gene expression (all conditions)') +
    ggtitle('Spaceflight data')
p
ggsave("plots/variance_mean_real_data_raw.pdf", plot=p, width=10, height=10, units='cm')

# If we were to use z-scores
# (Often done when using PCAs)
count_data_real_long_zscore <- count_data_real_long |>
  group_by(gene) |>
  mutate(counts_z = counts/sd(counts)) |>
  ungroup()
# Plot our gene
MYGENE = "AT2G41355"
p <- count_data_real_long_zscore |> 
    filter(gene==MYGENE) |>
    ggplot(aes(x=condition, y=counts_z))+
        geom_jitter(height = 0, width = 0.1) +
        theme_classic() +
        ggtitle(MYGENE)
p
# now calculate dispersions z
df_dispersion_real_zscore <- count_data_real_long_zscore |>
    group_by(gene) |>
    summarise(mean = mean(counts_z),
              sd   = sd(counts_z),
              variance = var(counts_z)) |>
    mutate(alpha = (variance-mean)/mean^2)
# now plot dispersions z
p <- ggplot(df_dispersion_real_zscore, aes(x=mean, y=variance)) +
    geom_point(color='grey') + 
    scale_x_log10() + scale_y_log10() + 
    theme_classic() +
    xlab('Mean gene expression (all conditions)') +
    ylab('Variance in gene expression (all conditions)') +
    ggtitle('Spaceflight data, z-scores')
p
ggsave("plots/variance_mean_real_data_zscore.pdf", plot=p, width=10, height=10, units='cm')
# similarly, we'd see a trend if we'd take the log
df_disperson_log <- count_data_real_long_zscore |>
    mutate(      
      log_counts   = log(counts+.1),
      log_counts_z = log(counts_z+.1)) |>
    group_by(gene) |>
    summarise(
      mean_z = mean(log_counts_z),
      sd_z   = sd(log_counts_z),
      variance_z = var(log_counts_z),
      mean = mean(log_counts), 
      sd = sd(log_counts),
      variance = var(log_counts))
# plot again
p <- ggplot(df_disperson_log, aes(x=mean, y=variance)) +
    geom_point(color='grey') + 
    scale_x_log10() + scale_y_log10() + 
    theme_classic() +
    xlab('Mean gene expression (all conditions)') +
    ylab('Variance in gene expression (all conditions)') +
    ggtitle('Spaceflight data')
p

################################################################################

# So now let's check out how dispersions look using vst
vsd <- vst(dds)
vst_counts <- assay(vsd)
df_dispersion_ddsvst <- vst_counts |>
    as.data.frame() |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "sample", values_to = "counts") |>
    group_by(gene) |>
    summarise(mean = mean(counts),
              sd   = sd(counts),
              variance = var(counts)) |>
    mutate(alpha = (variance - mean) / mean^2)
# now plot
p <- ggplot(df_dispersion_ddsvst, aes(x=mean, y=variance)) +
    geom_point(color='grey') + 
    scale_x_log10() + scale_y_log10() + 
    theme_classic() +
    xlab('Mean gene expression (all conditions)') +
    ylab('Variance in gene expression (all conditions)') +
    ggtitle('Spaceflight data, variance stabilized')
p
ggsave("plots/variance_mean_real_data_vst.pdf", plot=p, width=10, height=10, units='cm')

# now a pca plot by deseq2
pca_res <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
p <- ggplot(pca_res, aes(x=PC1, y=PC2, shape=condition,
                      color=condition)) +
    geom_point(size=3) +
    theme_classic() +
    xlab(paste0("PC1 (", round(100*pca_res$percentVar[1], 1), "%)")) +
    ylab(paste0("PC2 (", round(100*pca_res$percentVar[2], 1), "%)")) +
    ggtitle("PCA on vst data") +
    scale_color_manual(values=c("space_flight"="darkblue", "ground_control"="brown")) +
    theme(legend.position = "none")
p
ggsave("plots/pca_deseq2.pdf", plot=p, width=10, height=10, units='cm')

################################################################################
# Also in the overall DESeq2 analysis, dispersion-stablization is used
# (See also https://hbctraining.github.io/Intro-to-DGE/lessons/04b_DGE_DESeq2_analysis.html)

# So we have
# - Size factor correction
# - Dispersion correction

# Plot, remind ourselves that also the per-sample count normalization went well
normalized_counts_sampletotals <- colSums(normalized_counts)
p <- 
  data.frame(
    counts = unname(normalized_counts_sampletotals),
    sample = names(normalized_counts_sampletotals)
  ) |>
  mutate(sample = str_remove(sample, "Atha_WT_Col_0_sl_")) |>
  ggplot(aes(x = sample, y = counts)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Sample") + ylab("Counts\n(normalized by size factors)")
ggsave("plots/total_counts_per_sample_normalized.pdf", plot=p, width=7, height=10, units='cm')


# Run the whole analysis
dds_2 <- DESeq(dds) # saving to dds_2 here, usually not done
