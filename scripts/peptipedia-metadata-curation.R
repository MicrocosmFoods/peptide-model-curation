library(tidyverse)
library(UpSetR)

# cleaning and merging Peptipedia files
# read all individual peptidpedia CSV files into one df
peptipedia_dir <- "db_data/raw_data/peptipediadb"

# Get all CSV files that start with 'peptipedia-'
peptipedia_csvs <- list.files(peptipedia_dir, pattern=".csv", full.names = TRUE)

# Read and combine all CSV files
peptipedia_df <- peptipedia_csvs %>%
  map_df(~read_csv(.x) %>%
           mutate(bioactivity = str_remove(basename(.x), "peptipedia-") %>%
                    str_remove("\\.csv$")))
colnames(peptipedia_df) <- c("peptide_id", "sequence", "predicted", "bioactivity")

# pivot wider on bioactivity
peptipedia_df_wide <- peptipedia_df %>%
  select(-predicted) %>% 
  distinct(peptide_id, sequence, bioactivity) %>%
  mutate(value = TRUE) %>% 
  pivot_wider(names_from = bioactivity, values_from = value, values_fill = FALSE)

# pivot wider with just the validated sequences
peptipedia_validated_wide <- peptipedia_df %>% 
  filter(predicted == FALSE) %>% 
  select(-predicted) %>% 
  distinct(peptide_id, sequence, bioactivity) %>% 
  mutate(value = TRUE) %>% 
  pivot_wider(names_from = bioactivity, values_from = value, values_fill = FALSE)

# make upset plot of intersections of bioactivities
bioactivity_columns <- peptipedia_df_wide %>%
  select(-peptide_id, -sequence) %>%
  mutate(across(everything(), ~ case_when(
    . == TRUE ~ 1,
    . == FALSE ~ 0,
    TRUE ~ as.numeric(.)
  )))

validated_bioactivity_columns <- peptipedia_validated_wide %>% 
  select(-peptide_id, -sequence) %>%
  mutate(across(everything(), ~ case_when(
    . == TRUE ~ 1,
    . == FALSE ~ 0,
    TRUE ~ as.numeric(.)
  )))
  
# upset plot of overlap in bioactivities
upset(as.data.frame(bioactivity_columns), sets = c("anticancer", "antihypertensive", "antiinflammatory", "antimicrobial", "antioxidative", "immunomodulatory"), order.by = "freq")

upset(as.data.frame(validated_bioactivity_columns), sets = c("anticancer", "antihypertensive", "antiinflammatory", "antimicrobial", "antioxidative", "immunomodulatory"), order.by = "freq")

# read in peptipedia sequences from FermFooDB to see overlap with curated bioactivity
fermfoodb_peptipedia <- read.csv("db_data/raw_data/FermFooDB/peptipedia-fermfoodb-source.csv")
colnames(fermfoodb_peptipedia) <- c("peptide_id", "sequence")

# join with peptipedia wider df
fermfoodb_wide_bioactivities <- left_join(fermfoodb_peptipedia, peptipedia_df_wide)

# add column to peptipedia df that checks if record is in fermfoodb as well
peptipedia_df_wide <- peptipedia_df_wide %>%
  mutate(fermfoodb = peptide_id %in% fermfoodb_peptipedia$peptide_id)

peptipedia_validated_wide <- peptipedia_validated_wide %>% 
  mutate(fermfoodb = peptide_id %in% fermfoodb_peptipedia$peptide_id)

peptipedia_df_wide %>% 
  group_by(fermfoodb) %>% 
  count()

# write out TSVs

# metadata where a bioactivity can be experimentally/db validated or predicted with the Peptipedia models
write_tsv(peptipedia_df_wide, "db_data/cleaned_data/peptipedia/2024-11-04-peptipedia-metadata.tsv")

# metadata where TRUE values restricted to non-predicted (experimental evidence or in a DB)
write_tsv(peptipedia_validated_wide, "db_data/cleaned_data/peptipedia/2025-08-11-peptipedia-validated-metadata.tsv")

peptipedia_records <- peptipedia_df_wide %>% 
  select(peptide_id, sequence)

write_tsv(peptipedia_records, "db_data/cleaned_data/peptipedia/2024-11-04-peptipedia-records.tsv")
