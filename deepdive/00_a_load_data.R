

# %% ###########################################################################

# Loads two example datasets;
# - Synthetic one that was generated in main dir of this repo
# - Example dataset GLDS38, see https://osdr.nasa.gov/bio/repo/data/studies/OSD-38 (pre-processed to counts)

# this script yields
# - count_data_synth
# - count_data_real
# - count_data_synth_long
# - count_data_real_long


# %% ###########################################################################


# Import the table
count_data_synth <- read.table("data_frozen/count_table.tsv", 
                         sep = "\t", 
                         header = TRUE, 
                         row.names = 1, 
                         check.names = FALSE)
# Also import real data (spaceflight)
# (Downloaded from https://github.com/mishapaauw/RNA-seq-workshop/tree/main/data
# where Misha created a downsampled version of the original data.)
count_data_real <- read.csv("/Users/m.wehrens/Data_UVA/example-datasets/2026_RNAseq-Misha/GLDS38_raw_counts.csv",
                            row.names=1)
count_data_real_meta <- read.csv("/Users/m.wehrens/Data_UVA/example-datasets/2026_RNAseq-Misha/samples_to_condition.csv",
                          row.names=1)
print("Loaded count_data_synth and count_data_real")
                            

# make the format long
count_data_synth_long <- count_data_synth %>%
    rownames_to_column(var="gene") %>%
    pivot_longer(-gene, 
                names_to = "condition",
                values_to = "counts") 
    
count_data_real_long <- count_data_real %>%
    rownames_to_column(var="gene") %>%
    pivot_longer(-gene, 
                names_to = "sample",
                values_to = "counts") %>%
    # set condition to FLT or GC, based on presence of _FLT_ or _GC_ string in condtion column
    mutate(condition = case_when(
        str_detect(sample, "_FLT_") ~ "FLT",
        str_detect(sample, "_GC_") ~ "GC",
        TRUE ~ sample
    ))
print("Formatted data in long format, count_data_synth_long, count_data_real_long")
                