library(tidyverse)

# Modify output of feature count to match names of sample in summaryData.tsv
counts <- read_table("data/006/agc1_mm10.counts", skip = 1)

# Define function to compute geometric mean. If we add this number of counts to
# columns that don't exist we don't alter the geometric mean used by deseq2 to 
# compute size factors
geom_mean <- function(x){
  
  gm <- exp(mean(log(x)))
  return(gm)
  
}

counts_clean <- counts |> 
  select(-c(2:6)) |> 
  rename_with(.cols = 1:7, ~ c("ENSEMBL",
                               "wt_1",
                               "wt_2",
                               "wt_3",
                               "siagc1_1",
                               "siagc1_2",
                               "siagc1_3")) |> 
  
  rowwise() |> 
  mutate(wt_4 = ifelse(any(c(wt_1, wt_2, wt_3) == 0), 
                       0, 
                       round(geom_mean(c(wt_1, wt_2, wt_3)))),
         siagc1_4 = ifelse(any(c(siagc1_1, siagc1_2, siagc1_3) == 0), 
                           0, 
                           round(geom_mean(c(siagc1_1, siagc1_2, siagc1_3))))) |> 
  select(c("ENSEMBL",
            "wt_1",
            "wt_2",
            "wt_3",
            "wt_4",
            "siagc1_1",
            "siagc1_2",
            "siagc1_3",
            "siagc1_4"))

write_delim(counts_clean, "data/006/counts_table.tsv", delim = "\t")
