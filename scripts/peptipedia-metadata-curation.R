library(tidyverse)
library(UpSetR)

#################
# cleaning and merging Peptipedia files
#################

## First organize manually downloaded peptipedia files from the online API
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

## parse the PostgreSQL dumped and processed table of peptides with non-predicted, and non-other bioactivities
# read in and pivot wider so each peptide record is only present in the table once to make easier for FASTA generation and future metadata joining with results

peptipedia_postgresql_tsv <- read_tsv("db_data/raw_data/peptipedia_sql_dump/2025-11-25-non-predicted-bioactive-peptipedia.tsv")

peptipedia_postgresql_wide <- peptipedia_postgresql_tsv %>% 
  select(id_peptide, sequence, activity_name) %>% 
  mutate(value = TRUE) %>% # has this activity and is non-predicted moving forward, so the T/F gets inverted so the metadata joining makes sense in the context of comparing to BLAST-p results
  pivot_wider(
    names_from = activity_name, 
    values_from = value,
    values_fill = FALSE
  )

# write out this full wide metadata table
write_tsv(peptipedia_postgresql_wide, "db_data/cleaned_data/2025-11-25-cleaned-peptipedia-nonpredicted-metadata/2025-11-25-peptipedia-nonpredicted-all-bioactivities-records.tsv")

# output the first two columns for the peptide ID and sequence to conver to a FASTA
peptipedia_seq_records <- peptipedia_postgresql_wide %>% 
  select(id_peptide, sequence)

write_tsv(peptipedia_seq_records, "db_data/seqs/2025-11-25-non-predicted-records/2025-11-25-non-predicted-sequence-records.tsv")

# select bioactivities of interest 
peptipedia_metadata_filtered <- peptipedia_postgresql_wide %>% 
  select(id_peptide, sequence, `Angiotensin-converting enzyme (ace) inhibitors`, `Anti fungal`, `Anti inflamatory`, `Anti oxidative`, Antibacterial, Anticancer, Antimicrobial, Antiparasitic, Antiviral, `Blood brain barrier penetrating`, Cytotoxic, Immunological, Immunomodulatory, Immunoregulator, Metabolic, Neurological, Neuropeptide, Neurotoxin, `Signal peptide`, Toxic, Bacteriocin) %>% 
  mutate(peptide_id = id_peptide,
         anti_inflammatory = `Anti inflamatory`) %>% 
  select(peptide_id, sequence, everything()) %>% 
  select(-`Anti inflamatory`, -id_peptide)
  
# output select bioactivities metadata
write_tsv(peptipedia_metadata_filtered, "db_data/cleaned_data/2025-11-25-cleaned-peptipedia-nonpredicted-metadata/2025-11-25-peptipedia-nonpredicted-filtered-metadata.tsv")
