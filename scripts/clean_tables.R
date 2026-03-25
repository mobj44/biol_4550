library(tidyverse)
library(janitor)

df <- read_csv("data/tables/rodents_cytb_table.csv")
cleaned_df <- df %>%
    filter(!grepl("sp.", scientific_name)) %>%
    filter(!grepl("cf.", scientific_name)) %>%
    filter(!grepl("aff.", scientific_name)) %>%
    filter(!grepl("unclassified", scientific_name)) %>%
    filter(!grepl("mitochondrion", seq_name)) %>%
    filter(lengths(strsplit(as.character(scientific_name), " ")) > 1) %>%
    filter(seq_len > 400) %>%
    arrange(scientific_name)

by_len <- cleaned_df %>%
    arrange(desc(seq_len))
