#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 12.01.2021
# File last updated: 22.12.2022
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Tasks: Summarise transcript-level counts to gene level and create DESeq2
# data set for differential expression analysis.
#-------------------------------------------------------------------------------
# Inputs:
# Abundances from Kallisto pseudoalignment. Table describing the experimental
# design. Bter_v1_transcripts2genes.txt and virus_samples.csv files.

# Outputs:
# dds (DESeq2 data set) object.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(tximport)  # tximport()
library(DESeq2)  # DESeqDataSetFromTximport(), DESeq() 

# LOADING DATA ----
# Note:
# This function requires the "Bter_v1_transcripts2genes.txt" file and 
# "virus_samples.csv" file to be loaded in the environment.

# FUNCTION DEFINITION ----

BuildDDSForDEAnalysis <- function (x, y = "all_samples") {
  # Summarises transcript-level counts to gene level and create DESeq2
  # data set (dds).
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("brain", "fatbody" 
  #      or "ovaries").
  #   y: string denoting which set of samples should be analysed ("with_virus", 
  #      "no_virus" or "all_samples" (default)).
  #
  # Returns:
  #   A DESeq2 object with the transcript level counts summarised to gene level.
  
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
  
  # Set common internal variable.
  
  baseDir <- #"path/to/NER0008751_obj1_exp1_bter/02_outputs/"
  # This string is specific to your computer, PLEASE CHANGE accordingly.
  
  abundanceDir <- "22_kallisto_pseudoalignment_abundances"
  summaryDir <- "23_kallisto_pseudoalignment_summaries"
  study <- "_study_design.txt"
  
  # Generate character vector with sample names.
  
  sample_id <- 
    read.table(file.path(baseDir, tissuePath, summaryDir, 
                         paste0(x, study)),
               header = TRUE)
  
  # Add column containing full sample ID, to match the folder names.
  
  sample_id$full_id <- paste0(dir(file.path(baseDir, tissuePath, abundanceDir)))
  
  # Add condition column to sample_id.
  
  sample_id$condition <- paste0(sample_id$treatment, "_", sample_id$time_point) 
  
  # Add virus column to sample_id.
  
  sample_id <- merge(sample_id, virus_samples, by.x = "sample") 
  
  # Generate named vector with path to quantification files.
  
  files <- file.path(baseDir, tissuePath, abundanceDir, sample_id$full_id, "abundance.h5")
  
  names(files) <- paste0(sample_id$sample)
  
  # Import transcript-level estimates and summarise to gene level (Bter only)
  # producing "original counts and offsets".
  
  txi <- tximport(files, type = "kallisto", tx2gene = t2g)  

  # Initiate DESeq2 data set (dds) object.
  
  if (y == "all_samples") {
    ddsKallisto <- DESeqDataSetFromTximport(txi, sample_id, ~virus + condition)
    # Model has condition as the main effect, controlling for presence of 
    # virus-aligning samples.
  } else if (y == "with_virus" || y == "no_virus") {
    ddsKallisto <- DESeqDataSetFromTximport(txi, sample_id, ~condition)
    # Model has condition as the main effect.
  }
  
  # Collapse technical replicates.
  
  ddsKallisto <- collapseReplicates(ddsKallisto, ddsKallisto$sample, 
                                    ddsKallisto$full_id)
  
  # Remove fat body samples C_TP2G_rep4 and C_TP2G_rep6, which do not meet the
  # threshold of 12 million pseudoaligned reads with Kallisto.
  
  if (x == "fatbody") {
    ddsKallisto <- ddsKallisto[, -(c(10, 12))]
  }
  
  # Filter samples based on argument y.
  
  if (y == "no_virus") {
    ddsKallisto <- ddsKallisto[, -(which(ddsKallisto$virus == "Yes"))]
  } else if (y == "with_virus") {
    ddsKallisto <- ddsKallisto[, which(ddsKallisto$virus == "Yes")]
  }
  
  # Filter rows of dds object to remove any genes with less than 10 counts across 
  # all samples.
  
  ddsKallisto <- ddsKallisto[rowSums(counts(ddsKallisto)) >= 10, ]
  
  # Conduct differential expression analysis.
  
  ddsKallisto <- DESeq(ddsKallisto)
  
  # Return ddsKallisto
  
  return(ddsKallisto)
  
}
