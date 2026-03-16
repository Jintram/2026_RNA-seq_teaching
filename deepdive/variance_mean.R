
################################################################################

library(tidyverse)
library(ggplot2)

################################################################################

# Import the table
count_data <- read.table("data_frozen/count_table.tsv", 
                         sep = "\t", 
                         header = TRUE, 
                         row.names = 1, 
                         check.names = FALSE)

# make the format long
count_data_long <- count_data %>%
    rownames_to_column(var="gene") %>%
    pivot_longer(-gene, 
                names_to = "condition",
                values_to = "counts")
                
# now calculate mean, variance and sd
df_dispersion <- count_data_long %>%      
    group_by(gene) %>%
    summarise(mean = mean(counts),
              sd   = sd(counts),
              variance = var(counts)) %>%
    mutate(alpha = (variance-mean)/mean^2)

         

###

# DESeq2 has a dispersion estimate parameter, alpha, defined by
# var(K) = mu + alpha * mu^2
# Thus
# alpha = (var(K)-mu) / mu^2

### 
              
# Now plot variance versus mean
ggplot(df_dispersion, aes(x=mean, y=variance)) +
    geom_point()

# Now plot (var(K)-mu) / mu^2 versus mean
# This produces a similar plot to the dispresion estimate plot of DESeq2
ggplot(df_dispersion, aes(x=mean, y=alpha)) +
    geom_point() +
    scale_x_log10() + scale_y_log10(limits=c(1e-8, 10))+
    xlab("mean") + ylab("dispersion estimate\nalpha = (var(K)-mu) / mu^2")

###

# Note furthermore that DESeq2 selects the top-variance genes for taking 
# along in the PCA. 