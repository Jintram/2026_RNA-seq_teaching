
################################################################################

library(tidyverse)
library(ggplot2)

source("deepdive/00_a_load_data.R")
# - count_data_synth
# - count_data_real
# - count_data_synth_long
# - count_data_real_long
source("deepdive/00_b_calculate_dispersions.R")
# - df_dispersion_synth
# - df_dispersion_real

################################################################################
         

###

# DESeq2 has a dispersion estimate parameter, alpha, defined by
# var(K) = mu + alpha * mu^2
# Thus
# alpha = (var(K)-mu) / mu^2

### 
              
# Now plot variance versus mean
#
# Very nice source of information:
# https://hbctraining.github.io/DGE_workshop/lessons/01_DGE_setup_and_overview.html#rna-seq-count-distribution
# The red line gives the poisson relationship of the mean and variance
# The blue line gives the negative binominal relationship between mean and variance
fn_NBvariance <- function(mu, alpha) {
    return(mu + alpha * mu^2)
}
MY_ALPHA = 1/100
ggplot(df_dispersion, aes(x=mean, y=variance)) +
    geom_point() + 
    geom_abline(slope=1, intercept=0, color='red') +
    stat_function(fun=fn_NBvariance, args=list(alpha=MY_ALPHA), color='blue') +
    scale_x_log10() + scale_y_log10() + 
    theme_bw() +
    xlab('Mean gene expression (all conditions)') +
    ylab('Variance in gene expression (all conditions)')

# Now plot (var(K)-mu) / mu^2 versus mean
# This produces a similar plot to the dispresion estimate plot of DESeq2
ggplot(df_dispersion, aes(x=mean, y=alpha)) +
    geom_point() +
    scale_x_log10() + scale_y_log10(limits=c(1e-8, 10))+
    xlab("mean") + ylab("dispersion estimate\nalpha = (var(K)-mu) / mu^2")

###

# Note furthermore that DESeq2 selects the top-variance genes for taking 
# along in the PCA. 


### Now some plots for real data
# Let's first fit the dispersion coefficient
# Fit alpha in: variance = mu + alpha * mu^2
# Use nls on log-transformed values for a better fit across the wide range
df_fit <- df_dispersion_real %>% filter(mean > 0, variance > 0)
fit <- nls(log(variance) ~ log(mean + alpha * mean^2),
           data = df_fit,
           start = list(alpha = 0.01),
           control = nls.control(maxiter = 200))
MY_ALPHA_real <- coef(fit)[["alpha"]]
cat("Fitted alpha:", MY_ALPHA_real, "\n")
# 
p <- ggplot(df_dispersion_real, aes(x=mean, y=variance)) +
    geom_point(color='grey') + 
    geom_abline(slope=1, intercept=0, color='red', size=1) +
    stat_function(fun=fn_NBvariance, args=list(alpha=MY_ALPHA_real), color='blue', size=1) +
    scale_x_log10() + scale_y_log10() + 
    theme_classic() +
    annotate("text", x = 0.1, y = max(df_dispersion_real$variance)*.95, 
             label = paste0("variance = mean + alpha * mean^2\nalpha =",round(MY_ALPHA,2)),
             hjust = 0, vjust = 1) +
    # annotate("text", x = -Inf, y = Inf,
    #         label = "line1\nline2",
    #         hjust = 0, vjust = 1.5) +
    xlab('Mean gene expression (all conditions)') +
    ylab('Variance in gene expression (all conditions)') +
    ggtitle('Spaceflight data')
# p
# export to directory plots/
ggsave("plots/variance_mean_real_data.pdf", plot=p, width=10, height=10, units='cm')


