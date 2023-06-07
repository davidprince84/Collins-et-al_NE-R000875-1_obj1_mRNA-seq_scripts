#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Conduct gene ontology (GO) enrichment analysis on lists of 
# differentially expressed genes (DEGs) from the current study.
#-------------------------------------------------------------------------------
# Inputs:
# Overall gene lists for each tissue from the current study. OrthoFinder results.

# Outputs:
# .txt files stating the significantly enriched Biological Processes GO terms.
# NOTE: .txt files used rather than .csv files to prevent Excel from automatically
# formatting GeneRatio and BgRatio columns, which can result in dates being 
# inserted rather than fractions.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(clusterProfiler)  # enrichGO(), simplify().
library(org.Dm.eg.db)

# LOAD CUSTOM FUNCTIONS ----

source("f_IdentifySingleCopyOrthos.R")

# LOAD DATA ----

# Identify Bter and Dmel Single-Copy Orthologues ----

# Set working directory to OrthoFinder results.
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

setwd(paste0("../02_outputs/10_orthofinder/01_results/Orthologues/"))

# Identify single-copy orthologues between Dmel and Bter.

dmelBterOrthologues <- IdentifySingleCopyOrthos("Drosophila_melanogaster",
                                                "Bombus_terrestris")

# Rename columns.

colnames(dmelBterOrthologues) <- c("Orthogroup",
                                   "Dmel_gene_ID",
                                   "gene_symbol")

# FUNCTION DEFINITIONS ----

CalculateEnrichedGOTerms <- function(x, y, z) {
  # Calculates the enriched GO terms for a list of Bter DEGs using GO annotations
  # for the Dmel single-copy orthologues (as determined by OrthoFinder). Shows
  # which GO terms are redundant.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("brain", "fatbody" 
  #      or "ovaries").
  #   y: string denoting the name of treatment to be analysed. 
  #       ("C" or "R").
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in time point 2 compared to time 
  #      point one) or down-regulated genes (genes more highly expressed in time
  #      point 1 compared to time point 2) ("up" or "down").
  #
  # Returns:
  #   Nothing. A .txt file is saved to the working directory with the results.
  #
  # NOTE: The working directory needs to be set to the location where the results
  # are to be saved.
  
  # Set variables based on arguments.
  
  if (x == "brain") {
    allGeneColumnsToIndex <- c(1, 30, 31, 38, 39)
  } else if (x == "fatbody") {
    allGeneColumnsToIndex <- c(1, 30, 31, 38, 39)
  } else if (x == "ovaries") {
    allGeneColumnsToIndex <- c(1, 32, 33, 40, 41)
  } else {
    stop('Argument x must equal "brain", "fatbody" or "ovaries".')
  }
  
  if (y == "C") {
    allGenesWithOrthoscolumnToIndex <- c(1:3, 7)
    columnNameToIndex1 <- "DEG_in_C"
    columnNameToIndex2 <- "expression_with_age_in_C"
  } else if (y == "R") {
    allGenesWithOrthoscolumnToIndex <- c(1, 4, 5, 7)
    columnNameToIndex1 <- "DEG_in_R"
    columnNameToIndex2 <- "expression_with_age_in_R"
  } else {
    stop('Argument y must equal "C" or "R".')
  }
  
  if (z == "up") {
    directionOfDEGs <- "up-regulated"
  } else if (z == "down") {
    directionOfDEGs <- "down-regulated"
  } else {
    stop('Argument z must equal "up" or "down".')
  }
  
  # Load data.frame of all expressed genes for tissue x.
  
  setwd("../31_DESeq2_DEG_lists")
  
  allGenes <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  setwd("../40_GO_enrichment_analysis")
  
  # Remove columns no longer needed.
  
  allGenes <- allGenes[, allGeneColumnsToIndex]
  
  # Add Dmel Orthologues.
  
  allGenesWithOrthos <- merge(allGenes, dmelBterOrthologues)
  
  yGenesWithOrthos <- allGenesWithOrthos[, allGenesWithOrthoscolumnToIndex]
  
  # Index DEGs based on argument z.
  
  yDEGsWithOrthos <- 
    yGenesWithOrthos[yGenesWithOrthos[, columnNameToIndex1] == "YES" ,]
  
  yzDEGsWithOrthos <- 
    yDEGsWithOrthos[yDEGsWithOrthos[, columnNameToIndex2] == directionOfDEGs ,]
  
  # Conduct GO analysis.
  
  goEnrichResults <- enrichGO(yzDEGsWithOrthos$Dmel_gene_ID,
                              "org.Dm.eg.db",
                              keyType = "FLYBASE", 
                              ont = "BP",
                              universe = allGenesWithOrthos$Dmel_gene_ID,
                              qvalueCutoff = 0.05, minGSSize = 5)
  
  # Make data.frame of the enrichment results (if any).
  
  if (!is.null(goEnrichResults)) {
    goEnrichDF <- as.data.frame(goEnrichResults, stringsAsFactors = FALSE)
  }
  
  # Set file name.
  
  fileName <- paste(x, y, z, "regulated_DEGs_enriched_GO_results.txt", sep = "_") 
  
  # Save file with "No enriched GO terms found." if necessary, else
  # simplify GO results to remove redundancy.
  
  if (is.null(goEnrichResults)) {
    # Write a file stating that no enriched GO terms were returned.
    finalTable <- data.frame("Description" = "No enriched GO terms found.")
  } else if (length(goEnrichDF$ID) == 0) {
    # Write a file stating that no enriched GO terms were returned.
    finalTable <- data.frame("Description" = "No enriched GO terms found.")
  } else {
    # Simplify GO results.
    
    simpleGOResults <- as.data.frame(simplify(goEnrichResults),
                                     stringsAsFactors = FALSE)
    
    # Add column to specify terms are not redundant.
    
    simpleGOResults$redundant <- rep(0, times = length(simpleGOResults$ID))
    
    # Convert enrichResult object to data.frame.
    
    goEnrichDF <- as.data.frame(goEnrichResults, stringsAsFactors = FALSE)
    
    # Index the redundant terms from the results.
    
    redundantGODF <- goEnrichDF[!(goEnrichDF$ID %in% simpleGOResults$ID), ]
    
    # Add column specifying term is redundant.
    
    redundantGODF$redundant <- rep(1, times = length(redundantGODF$ID))
    
    # Bind the rows of the two data frames together.
    
    finalTable <- rbind(simpleGOResults, redundantGODF)
  } 
  
  # Save the results.
  
  write.table(finalTable, fileName, sep = "\t", row.names = FALSE)

} 

# EXECUTED STATEMENTS ----

# Calculate Enriched GO terms ----

# Brain.

# Set working directory.

setwd("../../../02_brain/")

# Create directory.

dir.create("40_GO_enrichment_analysis")
# Will produce a warning if directory already exists.

# Change directory.

setwd("40_GO_enrichment_analysis")

# Calculate enriched GO terms.

CalculateEnrichedGOTerms("brain", "C", "up")

CalculateEnrichedGOTerms("brain", "C", "down")

CalculateEnrichedGOTerms("brain", "R", "up")

CalculateEnrichedGOTerms("brain", "R", "down")

# Fat body.

# Create directory.

dir.create("../../03_fatbody/40_GO_enrichment_analysis")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../03_fatbody/40_GO_enrichment_analysis")

# Calculate enriched GO terms.

CalculateEnrichedGOTerms("fatbody", "C", "up")

CalculateEnrichedGOTerms("fatbody", "C", "down")

CalculateEnrichedGOTerms("fatbody", "R", "up")

CalculateEnrichedGOTerms("fatbody", "R", "down")

# Ovaries.

# Create directory.

dir.create("../../01_ovaries/40_GO_enrichment_analysis")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../01_ovaries/40_GO_enrichment_analysis")

# Calculate enriched GO terms.

CalculateEnrichedGOTerms("ovaries", "C", "up")

CalculateEnrichedGOTerms("ovaries", "C", "down")

CalculateEnrichedGOTerms("ovaries", "R", "up")

CalculateEnrichedGOTerms("ovaries", "R", "down")
