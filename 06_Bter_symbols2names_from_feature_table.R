#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Isolate gene symbols and corresponding gene names from features table.
#-------------------------------------------------------------------------------
# Inputs:
# Bombus terrestris (Bter) features table.

# Outputs:
# .csv file of B. terrestris gene symbols and corresponding gene names.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Feature Table ----

setwd("../00_data/05_feature_table")

featureTable <- read.delim("Bter_v1_feature_table.txt")

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Isolate Gene Symbols and Gene Names ----

# Index out the lines where feature = CDS.

allCDS <- featureTable[featureTable[, 1] == "CDS", ]

# Index out the lines where feature = ncRNA.

allncRNAs <- featureTable[featureTable[, 1] == "ncRNA", ]

# Index out the lines where feature = misc_RNA.

allmiscRNAs <- featureTable[featureTable[, 1] == "misc_RNA", ]

# Combine the CDS, ncRNA and misc_RNA data.fames together.

combinedFeatures <- rbind(allCDS, allncRNAs, allmiscRNAs)

# Remove duplication, so that a gene symbol is only attached to one name

noDuplicateSymbols <- combinedFeatures[!duplicated(combinedFeatures["symbol"]), ]
  
# Simplify data.frame to gene symbol and gene name only

genesToSymbols <- noDuplicateSymbols[, c("symbol", "name")]

# Format Gene Names ----

# Remove "isoform X etc" from gene names, as applies to transcripts, not genes.

for (ROW in 1:length(genesToSymbols$symbol)) {
  if (grepl("isoform X", genesToSymbols[ROW, 2])) {
    newName <- unlist(strsplit(genesToSymbols[ROW, 2], " isoform X"))
    genesToSymbols[ROW, 2] <- newName[1]
  }
}

# Remove "LOW QUALITY PROTEIN" from gene names.

for (ROW in 1:length(genesToSymbols$symbol)) {
  if (grepl("LOW QUALITY PROTEIN:", genesToSymbols[ROW, 2])) {
    newName <- unlist(strsplit(genesToSymbols[ROW, 2], "LOW QUALITY PROTEIN: "))
    genesToSymbols[ROW, 2] <- newName[2]
  }
}

# Write Data.Frame to File ----

# Save genesToSymbols to file.

write.csv(genesToSymbols, "genesymbols2genenames.txt",
          row.names = FALSE)
