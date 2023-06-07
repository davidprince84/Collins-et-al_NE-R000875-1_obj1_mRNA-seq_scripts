#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Generating indexes.
#-------------------------------------------------------------------------------
# Inputs:
# Bombus terrestris (Bter) GTF file and Bter and Holobee fasta files.
#
# Outputs:
# Indexes for Kallisto and HISAT2.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add HISAT2/2.1.0
module add kallisto/0.46.1

# STEP 1: UNZIP THE GENOME AND TRANSCRIPTOME FILES ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Change Directory ----

cd ../00_data/01_fasta

# Unzip Fasta Files ----

# Unzip genome fasta file.
# Args:
# --stdout: Keep file unchanged and write output to standard output.

gunzip --stdout Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz > Bter_v1_toplevel_dna.fasta

# Unzip transcriptome fasta file.

gunzip --stdout Bombus_terrestris.Bter_1.0.cdna.all.fa.gz > Bter_v1_cDNA.fasta

# STEP 2: UNZIP THE GTF FILE ----

# Change Directory ----

cd ../02_gtf

# Unzip GTF File ----

gunzip --stdout Bombus_terrestris.Bter_1.0.47.gtf.gz > Bter_v1.gtf

# STEP 3: GENERATE HISAT2 INDEX OF BTER GENOME ----

# Change Directory ----

cd ../01_fasta

# Build HISAT2 Index ----

hisat2-build Bter_v1_toplevel_dna.fasta Bter_v1_HISAT2_index

# STEP 4: GENERATE KALLISTO INDEX OF BTER TRANSCRIPTOME ----

# Indexing command.
# Args:
# -i: fileName for index to be created.

kallisto index -i Bter_v1_cDNA_index.idx Bter_v1_cDNA.fasta

# STEP 5: GENERATE KALLISTO INDEX OF HOLOBEE FASTA ----

kallisto index -i Holobee_v2016_1_index.idx HB_Bar_v2016.1.fasta
