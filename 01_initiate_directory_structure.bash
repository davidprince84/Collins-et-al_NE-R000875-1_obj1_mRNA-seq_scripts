#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Make a series of directories/subdirectories to conduct the analysis in.
# Move/copy files to directories.
#-------------------------------------------------------------------------------
# Inputs:
# The contents of the GitHub repository 
# https://github.com/davidprince84/Collins-et-al_NE-R000875-1_obj1_mRNA-seq_scripts.git
# The mRNA-seq read files from GEO GSE172422 downloaded using sra.
#
# Outputs:
# A series of named directories/subdirectories containing the scripts from
# the GitHub repository. Md5sums in the appropriate directory. 
# Gene lists in the appropriate directory.
#-------------------------------------------------------------------------------

# STEP 1: SET UP THE DIRECTORY STRUCTURE BY MAKING DIRECTORIES ----

# NOTE: This script assumes that the working directory is the downloaded GitHub repository.

# Making New Directories ----

# Make the initial directory for the analysis containing the project name.

mkdir ../NER0008751_obj1_bter

# Change into the new directory to make further directories.

cd ../NER0008751_obj1_bter

# Make directories to store the data, scripts and outputs.

mkdir 00_data
mkdir 01_scripts
mkdir 02_outputs

# Change into the data directory to make further directories.

cd 00_data

# Make directories to store different types of data.

mkdir 00_fastq
mkdir 01_fasta
mkdir 02_gtf
mkdir 03_bed
mkdir 04_sums
mkdir 05_feature_table
mkdir 06_gene_list_csv
mkdir 10_transcripts2genes
mkdir 11_virus_samples

# Change into the outputs directory to make further directories.

cd ../02_outputs

# Make directories to store outputs from analyses on all data, ovaries, head and fat body samples.

mkdir 00_all_tissues
mkdir 01_ovaries
mkdir 02_brain
mkdir 03_fatbody

# STEP 2: MOVE THE GITHUB REPOSITORY INTO THE SCRIPTS FOLDER ----

# Copy the Repository ----

# Change directory to the GitHub repository.

cd ../../Collins-et-al_NE-R000875-1_obj3_mRNA-seq_scripts

# Copy all the scripts and directories to the scripts folder.
# Args:
# -r: copy directories recursively

cp -r * ../NER0008751_obj1_bter/01_scripts

# Copy the .git directory to the scripts folder (as the previous command does not do this).

cp -r .git ../NER0008751_obj1_bter/01_scripts

# Delete the Repository ----

# Change directory.

cd ../

# Remove the empty repository. 
# Args:
# -r: remove directories and their contents recursively
# -f: do not ask permission to delete the files

rm -rf Collins-et-al_NE-R000875-1_obj3_mRNA-seq_scripts

# STEP 3: DELETE FILES AND SCRIPTS NOT NEEDED FOR THE MRNA-SEQ ANALYSIS ----

# Change Directory ----

cd NER0008751_obj1_bter

# Delete Files ----

rm -f 01_scripts/aggression.csv

rm -f 01_scripts/BORIS_video.csv

rm -f 01_scripts/cells.csv

rm -f 01_scripts/numbers.csv

rm -f 01_scripts/queen_measurements.csv

rm -f 01_scripts/survival.csv

rm -f 01_scripts/virus.csv

rm -f 01_scripts/NER0008751_Obj1_Exp1*.R  # Removes all R scripts not involved in the mRNA-seq analysis.

# STEP 4: MOVE MD5SUMS TO APPROPRIATE DIRECTORY ----

# Move Files ----

mv 01_scripts/md5sums*.txt 00_data/04_sums/

# STEP 5: MOVE GENE LISTS TO APPROPRIATE DIRECTORY ----

mv 01_scripts/pacifico_2018_table_s1_female.csv 00_data/06_gene_list_csv/

mv 01_scripts/chen_2014_table_s1.csv 00_data/06_gene_list_csv/

mv 01_scripts/korb_2021_table_s1.csv 00_data/06_gene_list_csv/

mv 01_scripts/kramer_2021_table_s3.csv 00_data/06_gene_list_csv/

# STEP 6: MOVE VIRUS SAMPLES DOCUMENT TO APPROPRIATE DIRECTORY ----

mv 01_scripts/virus_samples.csv 00_data/11_virus_samples/

# STEP 7: MOVE THE MRNA-SEQ READS ----

# Change directory to where the mRNA-seq reads were downloaded to.

# Use the following (commented out) code below to move the reads to the 
# 00_data/00_fastq folder where PATH/TO/ is the path from the current 
# directory to NER0008751_obj1_bter.

# mv *fastq.gz PATH/TO/NER0008751_obj1_bter/00_data/00_fastq  # PLEASE CHANGE.
