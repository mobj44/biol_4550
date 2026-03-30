library(tidyverse)
library(janitor)

cytb <- read_csv("data/tables/raw_muroidea_cytb_table.csv")
coi <- read_csv("data/tables/raw_muroidea_coi_table.csv")

cleaned_cytb <- cytb %>%
    filter(!grepl("sp.", scientific_name)) %>%
    filter(!grepl("cf.", scientific_name)) %>%
    filter(!grepl("aff.", scientific_name)) %>%
    filter(!grepl("unclassified", scientific_name)) %>%
    filter(lengths(strsplit(as.character(scientific_name), " ")) > 1) %>%
    arrange(scientific_name) %>%
    semi_join(coi, by="tax_id")


cleaned_coi <- coi %>%
    filter(!grepl("sp.", scientific_name)) %>%
    filter(!grepl("cf.", scientific_name)) %>%
    filter(!grepl("aff.", scientific_name)) %>%
    filter(!grepl("unclassified", scientific_name)) %>%
    filter(lengths(strsplit(as.character(scientific_name), " ")) > 1) %>%
    arrange(scientific_name) %>%
    semi_join(cytb, by="tax_id")

write_csv(cleaned_cytb, "data/tables/cleaned_muroidea_cytb_table.csv")
write_csv(cleaned_coi, "data/tables/cleaned_muroidea_coi_table.csv")


