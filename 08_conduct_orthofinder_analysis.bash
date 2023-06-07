#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Conduct OrthoFinder analysis for Bombus terrestris, Apis mellifera and
# Drosophila melanogaster.
#-------------------------------------------------------------------------------
# Inputs:
# B. terrestris, A. mellifera and D. melanogaster peptide sequences.  

# Outputs:
# OrthoFinder analysis results.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module add OrthoFinder/2.5.2

# STEP 1: CURATE DATA FOR ORTHOFINDER ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Copy the Relevant Data to a New Directory ----

# Change directory.

cd ../02_outputs

# Make new directory for the OrthoFinder analysis.

mkdir 10_orthofinder

# Copy peptide data for B. terrestris, A. mellifera and D. melanogaster 
# to the new directory.

cp ../00_data/01_fasta/*pep*.gz 10_orthofinder

# Change directory.

cd 10_orthofinder

# Unzip the files.

gunzip *.gz

# Loop over the peptide fasta files to extract the longest transcript per gene.
# Reason = quicker and more accurate (https://davidemms.github.io/orthofinder_tutorials/running-an-example-orthofinder-analysis.html)
# NOTE - path to primary_transcript.py is computer SPECIFIC, and therefore needs to be changed.

for FILE in *.fa
do
python #path/to/OrthoFinder-2.5.2/tools/primary_transcript.py $FILE  # CHANGE
done 

# Rename directory.

mv primary_transcripts 00_primary_transcripts

# Change directory.

cd 00_primary_transcripts

# Simplify file names so that they only refer to species.

mv Apis_mellifera.Amel_HAv3.1.pep.all.fa Apis_mellifera.fa
mv Bombus_terrestris.Bter_1.0.pep.all.fa Bombus_terrestris.fa
mv Drosophila_melanogaster.BDGP6.28.pep.all.fa Drosophila_melanogaster.fa

# Remove redudant files from 10_orthofinder directory.

# Change directory.

cd ../

# Remove files.

rm -f *pep*fa

# STEP 2: RUN ORTHOFINDER ----

# Run Orthofinder ----
# NOTE - path to orthofinder.py is computer SPECIFIC, and therefore needs to be changed.

# path/to/OrthoFinder-2.5.2/orthofinder.py -f 00_primary_transcripts  # CHANGE

# STEP 3: MOVE ORTHOFINDER RESULTS ----

# Rename Results Directory to Remove Date ----
# This simplifies later code, as there is then no need to change the date in the file path on 
# each run.

# Change directory.

cd 00_primary_transcripts/OrthoFinder/

# Rename directory.

mv Results* 01_results

# Move 01_results to 10_orthoFinder Directory ----
# Reduces length of directory paths in later code.

# Move 01_results directory.

mv 01_results ../../
