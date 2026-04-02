


################################################################################
# Show genes with padj < 0.05

res_ordered %>%
    as.data.frame() %>%
    filter(padj < 0.05)
    
    

################################################################################
# Make a plot of a specific gene

# Convert matrix to dataframe for plotting
df_plotting <- count_data %>% 
    rownames_to_column(var="gene") %>% 
    pivot_longer(- gene,
                 names_to = "condition",
                 values_to = "count") %>% 
    mutate(condition_strip = str_remove(condition, "\\.rep.*"))

df_plotting %>% 
    filter(gene=="aceA") %>% 
    # now plot
    ggplot(aes(x=condition_strip, y=count)) +
        geom_violin() +
        geom_point()
