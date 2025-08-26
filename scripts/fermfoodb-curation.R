library(tidyverse)

# cleaning and organizing FermFooDB peptide files
fermfoodb_path <- "db_data/raw_data/FermFooDB"
fermfoodb_files <- list.files(path = fermfoodb_path, pattern="\\.csv$", full.names = TRUE)

read_csv_files <- function(file) {
  data <- read.csv(file)
  substrate <- basename(file)
  data <- mutate(data, substrate = substrate)
  return(data)
}

all_fermfoodb <- do.call(rbind, lapply(fermfoodb_files, read_csv_files)) %>% 
  mutate(substrate = gsub("2024-fermfoodb-", "", substrate)) %>% 
  mutate(substrate = gsub("-peptides.csv", "", substrate))

all_fermfoodb %>% 
  group_by(Activity) %>% 
  count() %>% 
  arrange(desc(n)) %>% 
  print(n=50)

fermfoodb_seqs <- all_fermfoodb %>% 
  select(FMDB_ID, Sequence)

write_tsv(all_fermfoodb, "db_data/cleaned_data/FermFooDB/all-fermfoodb-peptide-metadata.tsv")
write_tsv(fermfoodb_seqs, "db_data/cleaned_data/FermFooDB/fermfoodb-seqs.tsv")
