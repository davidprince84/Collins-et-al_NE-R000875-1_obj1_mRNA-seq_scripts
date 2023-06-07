#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 10.12.2020
# File last updated: 21.12.2022
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Tasks: Summarise transcript-level counts to gene level and create DESeq2
# data set for quality control.
#-------------------------------------------------------------------------------
# Inputs:
# Abundances from Kallisto pseudoalignment. Table describing the experimental
# design.

# Outputs:
# dds (DESeq2 data set) object.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(tximport)  # tximport()
library(DESeq2)  # DESeqDataSetFromTximport(), DESeq() 

# LOADING DATA ----
# Note:
# This function requires the "Bter_v1_transcripts2genes.txt" file to be loaded 
# in the environment.

# FUNCTION DEFINITION ----

SumT2GBuildDDS <- function (x, y = "bter", z = "yes") {
  # Summarizes transcript-level counts to gene level and create DESeq2
  # data set.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("brain", "fatbody" 
  #      or "ovaries").
  #   y: string denoting whether the Kallisto abundances to be summarized are
  #      from the pseudoalignment to the Bombus terrestris (Bter) transcriptome
  #      or the Holobee sequences ("bter" (default) or "holobee").
  #   z: string denoting whether technical replicates should be collapsed 
  #      ("yes" (default) or "no").
  #
  # Returns:
  #   A DESeq2 object with the transcript level counts summarized to gene level.
  
  # Set variables based on argument.
  
  if (x == "brain") {
    tissuePath <- "02_brain"
  } else if (x == "fatbody") {
    tissuePath <- "03_fatbody"
  } else if (x == "ovaries") {
    tissuePath <- "01_ovaries"
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  if (y == "bter") {
    abundanceDir <- "22_kallisto_pseudoalignment_abundances"
    summaryDir <- "23_kallisto_pseudoalignment_summaries"
    study <- "_study_design.txt"
  } else if (y == "holobee") {
    abundanceDir <- "20_kallisto_pseudoalignment_holobee_abundances"
    summaryDir <- "21_kallisto_pseudoalignment_holobee_summaries"
    study <- "_holobee_study_design.txt"
  } else {
    stop('Argument y must be "bter" (default) or "holobee".')
  }
  
  # Set common internal variable.
  
  baseDir <- #"path/to/NER0008751_obj1_exp1_bter/02_outputs/"
  # This string is specific to your computer, PLEASE CHANGE accordingly.
  
  # Generate character vector with sample names ----
  
  sample_id <- 
    read.table(file.path(baseDir, tissuePath, summaryDir, 
                         paste0(x, study)),
               header = TRUE)
  
  # Add column containing full sample ID, to match the folder names.
  
  sample_id$full_id <- paste0(dir(file.path(baseDir, tissuePath, abundanceDir)))
  
  # Add condition column to sample_id.
  
  sample_id$condition <- paste0(sample_id$treatment, "_", sample_id$time_point) 
  
  # Generate named vector with path to quantification files.
  
  files <- file.path(baseDir, tissuePath, abundanceDir, sample_id$full_id, "abundance.h5")
  
  names(files) <- paste0(sample_id$sample)
  
  # Import transcript-level estimates and summarise to gene level (Bter only)
  # producing "original counts and offsets".
  
  if (y == "bter") {
    txi <- tximport(files, type = "kallisto", tx2gene = t2g)  
  } else if (y == "holobee") {
    txi <- tximport(files, type = "kallisto", txOut = TRUE)  
  } # As Holobee data can't be summarised to gene level.
  
  # Initiate DESeq2 data set (dds) object.
  
  ddsKallisto <- DESeqDataSetFromTximport(txi, sample_id, ~condition)
  
  # Filter rows of dds object to remove any genes with less than 10 counts across 
  # all samples.
  
  ddsKallisto <- ddsKallisto[rowSums(counts(ddsKallisto)) >= 10, ]
  
  # Collapse technical replicates if desired.
  
  if (z == "yes") {
    ddsKallisto <- collapseReplicates(ddsKallisto, ddsKallisto$sample, 
                                      ddsKallisto$full_id)
  } else if (z == "no") {
    ddsKallisto <- ddsKallisto
  } else {
    stop('Argument z must be "yes" (default) or "no".')
  }
  
  # Conduct differential expression analysis.
  
  ddsKallisto <- DESeq(ddsKallisto)
  
  # Return ddsKallisto
  
  return(ddsKallisto)
  
}
