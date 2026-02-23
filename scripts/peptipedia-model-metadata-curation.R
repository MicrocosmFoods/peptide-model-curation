library(tidyverse)

# curated metadata from non-predicted sequences for all activities
# from the raw SQL dump, created a table of all non-predicted sequences, their bioactivities, and the source of those sequences
peptipedia_non_predicted_sources <- read_tsv("db_data/raw_data/peptipedia_sql_dump/2025-11-25-non-predicted-bioactive-peptipedia.tsv")

# first get non-redundant set of the non-predicted sequences
# this is the peptipedia DB sequence file to compare to input sequences with DIAMOND Blast-p
nonpredicted_seqs <- peptipedia_non_predicted_sources %>% 
  select(id_peptide, sequence) %>% 
  mutate(peptide_id = id_peptide) %>%
  distinct(sequence, .keep_all = TRUE)

write_tsv(nonpredicted_seqs, "db_data/seqs/2026-02-11-nonpredicted-nonredundant-peptipedia-records.tsv")

# organize metadata
# keep only certain bioactivity categories for downstream comparison
bioactivity_list <- c("Angiotensin-converting enzyme (ace) inhibitors", "Anti fungal", "Anti oxidative", "Antibacterial", "Antimicrobial", "Antiparasitic", "Antiviral", "Blood brain barrier penetrating", "Cytotoxic", "Immunological", "Immunomodulatory", "Immunoregulator", "Metabolic", "Neurological", "Neuropeptide", "Neurotoxin", "Signal peptide", "Toxic", "Bacteriocin", "Anti inflamatory")

peptipedia_select_bioactivity_metadata <- peptipedia_non_predicted_sources %>% 
  select(id_peptide, activity_name) %>% 
  filter(activity_name %in% bioactivity_list) %>% 
  mutate(value = TRUE) %>% 
  pivot_wider(names_from = activity_name, values_from = value, values_fill = FALSE)

write_tsv(peptipedia_select_bioactivity_metadata, "db_data/cleaned_data/peptipedia/2025-11-25-cleaned-peptipedia-nonpredicted-metadata/2026-02-23-cleaned-peptipedia-nonpredicted-records-metadata.tsv")

# curatated antiinflammatory and immunomodulatory peptides from peptipedia for training datasets
# from the raw SQL dump, sourced antiinflammatory and immunomodulatory sequences that were non-predicted, have some inference from a training database or dataset, then reduced to distinct sequences/records
peptipedia_antiinflammatory <- read_tsv("db_data/raw_data/peptipedia_sql_dump/peptide_activity_sources_antiinflammatory.tsv")
peptipedia_immunomodulatory <- read_tsv("db_data/raw_data/peptipedia_sql_dump/peptide_activity_sources_immunomodulatory.tsv")

antiinflammatory_training_set <- peptipedia_antiinflammatory %>% 
  mutate(peptide_id = id_peptide) %>% 
  select(peptide_id, sequence) %>% 
  distinct(sequence, .keep_all = TRUE)

immunomodulatory_training_set <- peptipedia_immunomodulatory %>% 
  mutate(peptide_id = id_peptide) %>% 
  select(peptide_id, sequence) %>% 
  distinct(sequence, .keep_all = TRUE)


write_tsv(antiinflammatory_training_set, "db_data/model_training_data/2026-02-11-training-datasets/2026-02-11-antiinflammatory_training_set.tsv")
write_tsv(immunomodulatory_training_set, "db_data/model_training_data/2026-02-11-training-datasets/2026-02-11-immunomodulatory_training_set.tsv")

