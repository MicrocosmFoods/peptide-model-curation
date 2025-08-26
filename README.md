# Anti-inflammatory and Immunomodulatory Peptide Bioactivity Model and Metadata Curation

This repository contains raw data, scripts, and generated models for using machine learning classification models for predicting peptides with specific bioactivities. Publicly available databases were curated specifically for sequences and activities relevant to fermented foods, and models were generated for predicting peptides with anti-inflammatory and immunomodulatory activities. 

## Metadata curation

Peptides with associated metadata or curated bioactivities from the [Peptipedia] and [FermFooDB] databases were downloaded for further curation. The Peptipedia database does contain sequences from the FermFooDB database but doesn't retain important metadata features such as the food matrix the original peptide was obtained from. We downloaded each [food matrix CSV](https://webs.iiitd.edu.in/raghava/fermfoodb/food_matrix.php) from the FermFooDB and curated with the `scripts/fermfooddb-curation.R` script to organize verified bioactivities and other metadata. We then downloaded CSVs for select bioactivities from [Peptipedia](https://app.peptipedia.cl/download) where for each peptide is labelled with a column "Predicted" as True or False. We curated this metadata with `scripts/peptipedia-metadata-curation.R`. 

Since Peptipedia contains both labelled peptides and predicted peptides with their machine learning models for a select bioactivity, with curated sets of labelled seqeunces for generating machine learning classification models for anti-inflammatory and immunomodulatory activities with `scripts/peptipedia-model-curation.R`. We filtered to only sequences were Predicted = FALSE to obtain just labelled sequences that had some amount of validation or evidence for the bioactivity across the curated database sources collected through Peptipedia. Since our initial curation in Fall of 2024 Peptipedia no longer has an option to download sequences in TSV format and distinguishing in bulk whether a sequence has a labelled or predicted bioactivity. 

## Model generation

We generated models for predicting anti-inflammatory and immunomodulatory activities using the [AutoPeptideML v.1 webserver](http://peptide.ucd.ie/autopeptideml-1) and uploaded the TSVs in `db_data/model_training_data`. We selected "automatic search for negative peptides" and for building the anti-inflammatory model we selected immunomodulatory as an overlapping bioactivity and vice versa for the immunomdulatory model. We kept the homology partitioning threshold at the default of 0.3 and using the default alignment algorithm as mmseqs. Model output files are located in `models/` and split by bioactivity. 

For the anti-inflammatory model we created two models, ANIF_1 and ANIF_2. ANIF_1 was trained using labelled sequences from the Peptipedia database as described above. ANIF_2 was constructed using positive and negative labelled sequences from a benchmark database collected Immune Epitope Database curated as part of the AIPpred and AntiInflam machine learning models. 

## Datasets

We predicted the anti-inflammatory and immunomodulatory potential of peptides predicted from genomes assembled from fermented foods and proteomics experiments from fermented foods, as well as predicted other bioactivities with models from the AutoPeptideML release. The datasets as well as machine learning models and benchmarking sequences are also available [on Zenodo](https://zenodo.org/records/16749254).