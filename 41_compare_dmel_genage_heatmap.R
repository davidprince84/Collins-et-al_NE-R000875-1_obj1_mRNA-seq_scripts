#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Comparison of Bombus terrestris (Bter) differentially expressed genes 
# (DEGs) with Drosophila melanogaster (Dmel) GenAge data set.
#-------------------------------------------------------------------------------
# Inputs:
# GenAge Dmel gene list with Bter orthologs (from script #40). Lists of 
# Bter DEGs.

# Outputs:
# A heatmap of the relevant results.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# org.Dm.eg.db loaded via the script in the "Load data" section.
# pheatmap and ggplot2 loaded via the script in the "Custom functions" section.

# LOAD DATA ----

# Load GenAge Data ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("40_compare_dmel_genage_prep.R")

# Add a new column for Bter expression to genAgeOrtho.

genAgeOrtho$Study_expression <- NA

# LOAD CUSTOM FUNCTIONS ----

source("f_PlotCustomHeatmap.R")

# FUNCTION DEFINITIONS ----

ReportPresenceOfGenAgeGenes <- function (x, y = "none", z = "none") {
  # Reports whether a list of DEGs contains orthologues of the Dmel GenAge genes.
  #
  # Args:
  #   x: string denoting the study that the gene list is taken from ("current", 
  #      "chen" or "pacifico").
  #   y: string denoting which tissue is to be analysed for Bter genes from the
  #      current study ("none" (default), "brain", "fatbody" or "ovaries").
  #   z: string denoting the name of treatment from the current study to be analysed. 
  #      ("none" (default), "C" or "R").
  #
  # Returns:
  #   genAgeOrtho data.frame with presence of differentially expressed genes 
  #   for stated study added.
  
  # Load DEGs list based on arguments.
  
  if (x == "current") {
    upregList <- read.csv(paste0(y, "_results_", z, 
                                 "_treatment_up_regulated_LFC0.csv"))
    downregList <- read.csv(paste0(y, "_results_", z, 
                                   "_treatment_down_regulated_LFC0.csv"))  
  } else if (x == "chen") {
    upregList <- chenUpDEGsFB
    downregList <- chenDownDEGsFB
  } else if (x == "pacifico") {
    upregList <- pacificoUpDEGsFB
    downregList <- pacificoDownDEGsFB
  } else {
    stop('Argument x must be "current", "chen" or "pacifico".')
  }
  
  # Designate internal copy of genAgeOrtho.
  
  genAge <- genAgeOrtho
  
  # Set variables based on arguments.
  
  if (x == "current") {
    speciesGeneID <- "Bter_gene_ID"
    listGeneID <- "symbol"
  } else {
    speciesGeneID <- "Flybase_gene_ID"
    listGeneID <- "Dmel_gene_ID"
  }
  
  # Loop over GenAge data and add "1" to Study_expression column.
  
  for(i in (1:length(genAge[, 1]))) {
    if (genAge[i, speciesGeneID] %in% upregList[, listGeneID]) {
      genAge[i, "Study_expression"] <- 1
    } else if (genAge[i, speciesGeneID] %in% downregList[, listGeneID]) {
      genAge[i, "Study_expression"] <- -1
    } else {
      genAge[i, "Study_expression"] <- 0
    }
  }
  
  # Return genAge.
  
  return(genAge)
  
}

# EXECUTED STATEMENTS ----

# Bter Current Study ----

# Ovaries.

# Change working directory.

setwd("../02_outputs/01_ovaries/31_DESeq2_DEG_lists/")

# Run function.

ovariesControlResults <- ReportPresenceOfGenAgeGenes("current", "ovaries", "C")

ovariesRemovalResults <- ReportPresenceOfGenAgeGenes("current", "ovaries", "R")

# Fat body.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

fatBodyControlResults <- ReportPresenceOfGenAgeGenes("current", "fatbody", "C")

fatBodyRemovalResults <- ReportPresenceOfGenAgeGenes("current", "fatbody", "R")

# Brain.

setwd("../../02_brain/31_DESeq2_DEG_lists/")

brainControlResults <- ReportPresenceOfGenAgeGenes("current", "brain", "C")

brainRemovalResults <- ReportPresenceOfGenAgeGenes("current", "brain", "R")

# Dmel Previous Studies ----

# Chen et al. (2014).

chenResults <- ReportPresenceOfGenAgeGenes("chen")

# Pacifico et al. (2018).

pacificoResults <- ReportPresenceOfGenAgeGenes("pacifico")

# Combine results.

combinedResults <- genAgeOrtho[, c(1, 3, 4, 9, 11)]

# Add current and previous study results.

combinedResults$pacifico_dmel_brain <- pacificoResults[, 12]

combinedResults$current_brain_removal <- brainRemovalResults[, 12]

combinedResults$current_brain_control <- brainControlResults[, 12]

combinedResults$chen_dmel_fatbody <- chenResults[, 12]

combinedResults$current_fatbody_removal <- fatBodyRemovalResults[, 12]

combinedResults$current_fatbody_control <- fatBodyControlResults[, 12]

combinedResults$current_ovaries_removal <- ovariesRemovalResults[, 12]

combinedResults$current_ovaries_control <- ovariesControlResults[, 12]

# Plot and Save Heatmap ----

# Prepare data for heatmap.

# Heatmap showing genes which are differentially expressed in a Dmel study 
# (Chen or Pacifico) and at least 1 Bter sample, or in Bter brain or fat body.

# Remove rows where only DEGs are in ovaries.

noBterZeroResults <- combinedResults[rowSums(combinedResults[, c(7, 8, 10, 11)] == 0, na.rm=TRUE) < 4, ]

# Remove ‘Juvenile hormone esterase binding protein 29 (DmP29)’ as it is annotated
# with both pro- and anti-longevity influences, and R won't allow
# this to be annotated on a heatmap (as the row names in the annotation
# data.frame would be the same).

noBterZeroResults <- noBterZeroResults[!(noBterZeroResults[, "name"] == "Juvenile hormone esterase binding protein 29 (DmP29)"), ]

# Format names where necessary.

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0003028"), "name"] <- "ovo" 

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0020238"), "name"] <- "14-3-3 epsilon" 

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0023215"), "name"] <- "Mnt" 

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0023511"), "name"] <- "ER degradation enhancer, mannosidase alpha-like 1"

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0026317"), "name"] <- "Tuberous sclerosis complex genes 1"

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0030512"), "name"] <- "NAD synthetase"

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0036813"), "name"] <- "Autophagy-related 3"

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0033153"), "name"] <- "Growth arrest and DNA damage-inducible 45"

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0260990"), "name"] <- "yata"

noBterZeroResults[which(noBterZeroResults$Flybase_gene_ID == "FBgn0262169"), "name"] <- "magu"

# Plotting matrix of +1/0/-1 values.

plotMatrix <- as.matrix (noBterZeroResults[, c(6:13)])

# Add gene names to the matrix.

row.names(plotMatrix) <- noBterZeroResults[, "name"]

# Plot heatmap.

heatmapPlot <- PlotCustomHeatmap(41)

# Save heatmap.

setwd("../../00_all_tissues")

ggsave("30_NER0008751_obj1_fig_S14_GenAge_heatmap.svg", plot = heatmapPlot, 
       width = 21, height = 29, units = "cm")
