#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Download GenAge model organism data set and select Drosophila
# melanogaster genes from it.
#-------------------------------------------------------------------------------
# Inputs:
# None.

# Outputs:
# .csv file containing only the Drosophila melanogaster genes from the 
# GenAge database.
#-------------------------------------------------------------------------------

# STEP 1: DOWNLOAD GENAGE FILE ----

# Change Directory ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

cd ../00_data/06_gene_list_csv

# Download Model Organisms Data From GenAge Website ----

wget https://genomics.senescence.info/genes/models_genes.zip

# Unzip the Downloaded Folder and Remove it ----

unzip models_genes.zip && rm -f models_genes.zip

# Rename Release Notes to Make More Informative ----

mv release.html genage_release.html 

# STEP 2: SELECT DROSOPHILA MELANOGASTER GENES ----

# Save Lines with Drosophila melanogaster Genes to a new File ----

grep "GenAge ID" genage_models.csv > Dmel_genage.csv  # save the headers.

grep "Drosophila melanogaster" genage_models.csv >> Dmel_genage.csv  # add the genes.
