library(tidyverse)
library(janitor)

cytb <- read_csv("data/tables/rodents_cytb_table.csv")
coi <- read_csv("data/tables/rodents_coi_table.csv")

cleaned_cytb <- cytb %>%
    filter(!grepl("sp.", scientific_name)) %>%
    filter(!grepl("cf.", scientific_name)) %>%
    filter(!grepl("aff.", scientific_name)) %>%
    filter(!grepl("unclassified", scientific_name)) %>%
    filter(!grepl("mitochondrion", seq_name)) %>%
    filter(lengths(strsplit(as.character(scientific_name), " ")) > 1) %>%
    filter(seq_len > 400) %>%
    arrange(scientific_name) %>%
    semi_join(coi, by="tax_id")


cleaned_coi <- coi %>%
    filter(!grepl("sp.", scientific_name)) %>%
    filter(!grepl("cf.", scientific_name)) %>%
    filter(!grepl("aff.", scientific_name)) %>%
    filter(!grepl("unclassified", scientific_name)) %>%
    filter(!grepl("mitochondrion", seq_name)) %>%
    filter(lengths(strsplit(as.character(scientific_name), " ")) > 1) %>%
    filter(seq_len > 400) %>%
    arrange(scientific_name) %>%
    semi_join(cytb, by="tax_id")

write_csv(cleaned_cytb, "data/tables/cleaned_rodents_cytb_table.csv")
write_csv(cleaned_coi, "data/tables/cleaned_rodents_coi_table.csv")
