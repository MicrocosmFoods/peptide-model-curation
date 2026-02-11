library(tidyverse)

# curating antiinflammatory and immunomodulatory peptides from peptipedia for training datasets
# from the raw SQL dump, sourced antiinflammatory and immunomodulatory sequences that were non-predicted, have some inference from a training database or dataset
peptipedia_antiinflammatory <- read_tsv("db_data/raw_data/peptipedia_sql_dump/peptide_activity_sources_antiinflammatory.tsv")
peptipedia_immunomodulatory <- read_tsv("db_data/raw_data/peptipedia_sql_dump/peptide_activity_sources_immunomodulatory.tsv")

antiinflammatory_training_set <- peptipedia_antiinflammatory %>% 
  mutate(peptide_id = id_peptide) %>% 
  select(peptide_id, sequence) %>% 
  distinct(sequence, .keep_all = TRUE)

immunomodulatory_training_set <- peptipedia_immunomodulatory %>% 
  mutate(peptide_id = id_peptide) %>% 
  select(peptide_id, sequence)


write_tsv(antiinflammatory_training_set, "db_data/model_training_data/antiinflammatory_training_set.tsv")
write_tsv(immunomodulatory_training_set, "db_data/model_training_data/immunomodulatory_training_set.tsv")

