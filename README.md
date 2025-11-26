# Peptipedia Database Metadata Curation and Custom Model Creation

This repository contains documentation, raw data, scripts, and generated models for using machine learning classification models for predicting peptides with specific bioactivities. Publicly available databases were curated specifically for sequences and activities relevant to fermented foods, and models were generated for predicting peptides with anti-inflammatory and immunomodulatory activities. 

## Important Note

When Peptipedia was first released and published sometime in 2024, the raw database files were not available. Files that were available to be manually downloaded from the Peptipedia database website had to be done so by individual activity or original database source. Originally Peptipedia allowed downloads of peptides with a specific activity to be downloaded in FASTA or TSV format, importantly where the TSV format contained a "predicted" column with `True` or `False` entries. This was an important distinction, as sequences with the `False` labels for prediction are those that have some sort of evidence from another curated database or experimental evidence in the literature of that bioactivity. Whereas peptides with the `True` labels for that bioactivity were predicted using models developed by Peptipedia. It's important to note that even sequences with a `False` predicted label could still be technically predicted by a different tool/database, but were not done so by Peptipedia's models as part of their publication, which predicts a huge amount of sequences. Including these sequences in sequence-based comparisons and especially devloping custom models below could lead to false-positives or worse data leakage and overfitting for generated models. However, because this could only be done manually at the time and we were only interested in specific bioactivities, we downloaded only a select few of these TSV files which are listed in [db_data/raw_data/peptipediadb/](db_data/raw_data/peptipediadb/), all ending in .csv. and have the structure of `id_peptide, sequence, predicted`.

However, at some point in 2025 the online Peptipedia database no longer had an option to download the TSV file for a given peptide bioactivity, and one could no longer filter sequences that were not predictions. However, the raw database file was made available by the original authors during this time and is available to download [through Google Drive](https://drive.google.com/file/d/1uvTGOdjpsPYxvx00g8KbMv5tTDKsjSAg/view?usp=drive_link). This file is the PostgreSQL dump of the Peptipedia database. The next section documents how the database was downloaded, processed, and transformed to obtain a necessary TSV file containing the "predicted" labels for all peptide sequences regardless of bioactivity.

## Accessing and Processing Raw PostgreSQL Peptipedia Database Dump

### Installation and Setup

First, install PostgreSQL if you don't already have it on your system. I'm on a Macbook Pro with a Apple M3 Max Chip, so I installed PostgreSQL16 with homebrew and added it to my path with:

```
arch -arm64 brew install postgresql@16
echo 'export PATH="/opt/homebrew/opt/postgresql@16/bin:$PATH"' >> ~/.zshrc
```

Then start the server and ensure proper installation with:
```
brew services start postgresql@16
psql --version
```

If the last command returns a version, you should be ready for the next step, which is to create a database and load the PostgreSQL dumped file into it with:

```
createdb peptides_db
psql -d peptides_db -f raw_data/peptipedia.20240323.sql
```

Then start the server with:

```
psql peptides_db
```

You should then be able to see all the loaded datatables by entering `\dt`: 

```
peptides_db=# \dt
                       List of relations
 Schema |               Name               | Type  
--------+----------------------------------+-------
 public | activity                         | table
 public | gene_ontology                    | table 
 public | peptide                          | table
 public | peptide_has_activity             | table
 public | peptide_has_go                   | table
 public | peptide_has_pfam                 | table
 public | peptide_has_source               | table
 public | pfam                             | table
 public | predictive_model                 | table
 public | source                           | table
(13 rows)
```

Now we are ready to process the data tables so it's usable for our downstream purposes of joining peptide sequences for only non-predicted activities, which is the necessary information for us to make DIAMOND-blastp comparisons to our genomic sequences of interest and create training datasets for machine learning models. 

### Create a new table with activity labels and sequences
The information for the peptide activity, whether the activity is predicted, and the sequence of the peptide is contained in three different tables. The commands below create subset datatables to join that information, and then filter to just sequences where predicted = False.

First join peptide activity labels to the peptide list:
```
CREATE TABLE peptide_activity_labels AS
SELECT
  peptide_has_activity.id_peptide,
  peptide_has_activity.id_activity,
  activity.name AS activity_name,
  peptide_has_activity.predicted
FROM peptide_has_activity
JOIN activity
  ON peptide_has_activity.id_activity = activity.id_activity;
```

Then join with the peptide sequences:
```
DROP TABLE IF EXISTS peptide_activity_labels_with_seq;

CREATE TABLE peptide_activity_labels_with_seq AS
SELECT
  peptide_activity_labels.*,
  peptide.sequence
FROM peptide_activity_labels
JOIN peptide
  ON peptide.id_peptide = peptide_activity_labels.id_peptide;
```

Then filter down to only rows where the prediction column = False, and write out to a TSV file:
```
DROP TABLE IF EXISTS peptide_activity_not_predicted;

CREATE TABLE peptide_activity_not_predicted AS
SELECT *
FROM peptide_activity_labels_with_seq  -- or your current table name
WHERE predicted = false;

\copy peptide_activity_not_predicted TO 'peptide_activity_not_predicted.tsv' WITH (FORMAT csv, HEADER true, DELIMITER E'\t');
```

This TSV was then used as the basis for creating a main metadata file for non-predicted (as least not by Peptipedia directly) peptide activities and a FASTA file for DIAMOND-blastp sequence comparisons. Additionally the peptide_id matches the Peptipedia record ID that can be searched through their web portal to browse what other databases the peptide was sourced from and any cited literature the peptide and that activity may have been originally sourced from as well. Following up on this information is important when assessing the validity of models and DIAMOND-Blastp sequence-based comparison hits, and you should now be able to easily find that information with the metadata parsed in this fashion.

The cleaned peptipedia metadata for peptides that are labelled as non-predicted and non-other is in `db_data/cleaned_data/2025-11-25-cleaned-peptipedia-nonpredicted-metadata` with the complete TSV of the transformed metadata where each peptide record is one row, and each column is TRUE or FALSE for that bioactivity. In this table note that the TRUE/FALSE is inverted from before. In this case it's TRUE if it has that bioactivity, and FALSE if not. This file is to match up with results such as DIAMOND-BLASTp results and therefore make the TRUE/FALSE case more intuitive for viewing results like that. The FASTA file for these sequences are in `db_data/seqs/2025-11-25-non-predicted-records`. 

## 2024 Metadata Curation

Peptides with associated metadata or curated bioactivities from the Peptipedia and FermFooDB databases were downloaded for further curation. The Peptipedia database does contain sequences from the FermFooDB database but doesn't retain important metadata features such as the food matrix the original peptide was obtained from. We downloaded each [food matrix CSV](https://webs.iiitd.edu.in/raghava/fermfoodb/food_matrix.php) from the FermFooDB and curated with the `scripts/fermfooddb-curation.R` script to organize verified bioactivities and other metadata. We then downloaded CSVs for select bioactivities from [Peptipedia](https://app.peptipedia.cl/download) where for each peptide is labelled with a column "Predicted" as True or False. We curated this metadata with `scripts/peptipedia-metadata-curation.R`. 

Since Peptipedia contains both labelled peptides and predicted peptides with their machine learning models for a select bioactivity, with curated sets of labelled seqeunces for generating machine learning classification models for anti-inflammatory and immunomodulatory activities with `scripts/peptipedia-model-curation.R`. We filtered to only sequences were Predicted = FALSE to obtain just labelled sequences that had some amount of validation or evidence for the bioactivity across the curated database sources collected through Peptipedia. Since our initial curation in Fall of 2024 Peptipedia no longer has an option to download sequences in TSV format and distinguishing in bulk whether a sequence has a labelled or predicted bioactivity. 

## Model Generation for Anti-inflammatory and Immunomodulatory Activities

We generated models for predicting anti-inflammatory and immunomodulatory activities using the [AutoPeptideML v.1 webserver](http://peptide.ucd.ie/autopeptideml-1) and uploaded the TSVs in `db_data/model_training_data`. We selected "automatic search for negative peptides" and for building the anti-inflammatory model we selected immunomodulatory as an overlapping bioactivity and vice versa for the immunomdulatory model. We kept the homology partitioning threshold at the default of 0.3 and using the default alignment algorithm as mmseqs. Model output files are located in `models/` and split by bioactivity. 

For the anti-inflammatory model we created two models, ANIF_1 and ANIF_2. ANIF_1 was trained using labelled sequences from the Peptipedia database as described above. ANIF_2 was constructed using positive and negative labelled sequences from a benchmark database collected Immune Epitope Database curated as part of the AIPpred and AntiInflam machine learning models. 

## Datasets

We predicted the anti-inflammatory and immunomodulatory potential of peptides predicted from genomes assembled from fermented foods and proteomics experiments from fermented foods, as well as predicted other bioactivities with models from the AutoPeptideML release. The datasets as well as machine learning models and benchmarking sequences are also available [on Zenodo](https://zenodo.org/records/16749254).