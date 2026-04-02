
#%% ###############################################################################

source("deepdive/00_a_load_data.R")
# - count_data_synth
# - count_data_real
# - count_data_synth_long
# - count_data_real_long

library(tidyverse)
library(ggplot2)

library(DESeq2)

#%% ################################################################################
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


#%% ################################################################################
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


#%% ################################################################################
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
p
ggsave("plots/pca_deseq2.pdf", plot=p, width=10, height=10, units='cm')

#%% ################################################################################
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
