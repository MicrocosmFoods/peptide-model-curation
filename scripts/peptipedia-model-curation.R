library(tidyverse)

# curating antiinflammatory, immunomodulatory, and antioxidative seqs
# selecting predicted = FALSE for training models
peptipedia_antiinflammatory <- read.csv("db_data/raw_data/peptipediadb/peptipedia-antiinflammatory.csv")
peptipedia_immunomodulatory <- read.csv("db_data/raw_data/peptipediadb/peptipedia-immunomodulatory.csv")

antiinflammatory_training_set <- peptipedia_antiinflammatory %>% 
  filter(predicted == 'False') %>% 
  mutate(peptide_id = id_peptide) %>% 
  select(peptide_id, sequence)

immunomodulatory_training_set <- peptipedia_immunomodulatory %>% 
  filter(predicted == 'False') %>% 
  mutate(peptide_id = id_peptide) %>% 
  select(peptide_id, sequence)


write_tsv(antiinflammatory_training_set, "db_data/model_training_data/antiinflammatory_training_set.tsv")
write_tsv(immunomodulatory_training_set, "db_data/model_training_data/immunomodulatory_training_set.tsv")

