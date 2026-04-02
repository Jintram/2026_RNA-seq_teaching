

# this script yields
# - df_dispersion_synth
# - df_dispersion_real

source("deepdive/00_a_load_data.R")
# - count_data_synth
# - count_data_real
# - count_data_synth_long
# - count_data_real_long


# now calculate mean, variance and sd
df_dispersion_synth <- count_data_synth_long %>%      
    group_by(gene) %>%
    summarise(mean = mean(counts),
              sd   = sd(counts),
              variance = var(counts)) %>%
    mutate(alpha = (variance-mean)/mean^2)
df_dispersion_real <- count_data_real_long %>%      
    group_by(gene) %>%
    summarise(mean = mean(counts),
              sd   = sd(counts),
              variance = var(counts)) %>%
    mutate(alpha = (variance-mean)/mean^2)
print("Calculated dispersion data frame, df_dispersion_synth, df_dispersion_real")
         