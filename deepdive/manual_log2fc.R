

source("deepdive/import_data.R")

library(tidyverse)
library(ggplot2)

# Calculate means and log2fc per condition
count_data_real_long_l2fc <- count_data_real_long |>
    # first gather averages
    group_by(condition, gene) |>
    summarise(gene_average = mean(counts)) |>    
    # now determine FC
    filter(gene_average>0) |>
    group_by(gene) |>
    pivot_wider(names_from = condition, values_from = gene_average) |>
    mutate(FC = FLT/GC, 
            L2FC = log2(FC)) |>
    arrange(desc(L2FC))
    


# now make a small plot for gene="AT2G41355"
MYGENE = "AT2G41355"
p <- count_data_real_long |> 
    filter(gene==MYGENE) |>
    ggplot(aes(x=condition, y=counts))+
        geom_jitter(height = 0, width = 0.1) +
        theme_classic() +
        ggtitle(MYGENE)
ggsave("plots/gene_example_highl2fc.pdf", plot=p, width=7, height=7, units='cm')


