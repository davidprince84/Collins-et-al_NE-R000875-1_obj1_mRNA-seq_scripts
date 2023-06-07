#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Compare Bombus terrestris (Bter) differentially 
# expressed genes (DEGs) with Drosophila melanogaster (Dmel) enzymatic antioxidant
# network genes, as defined by Kramer et al. (2021).
#-------------------------------------------------------------------------------
# Inputs:
# Bter DEGs from the current study. OrthoFinder results. List of Dmel enzymatic 
# antioxidant network genes from Kramer et al. (2021).

# Outputs:
# Heatmap of enzymatic antioxidant network genes differentially expressed in the
# Bter DEGs.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# pheatmap and ggplot2 loaded via the script in the "Custom functions" section.

# LOAD CUSTOM FUNCTIONS ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("f_IdentifySingleCopyOrthos.R")

source("f_PlotCustomHeatmap.R")

# LOAD DATA ----

# Load Kramer et al. (2021) Data ----

# Set working directory.

setwd("../00_data/06_gene_list_csv")

# Load data.

antioxidantGenes <- read.csv("kramer_2021_table_s3.csv",
                           stringsAsFactors = FALSE)

# Set working directory to OrthoFinder results.

setwd("../../02_outputs/10_orthofinder/01_results/Orthologues/")

# Identify single-copy orthologues between Dmel and Bter.

dmelBterOrthologues <- IdentifySingleCopyOrthos("Drosophila_melanogaster",
                                                "Bombus_terrestris")

# Rename columns of OrthoFinder results to allow merge.

colnames(dmelBterOrthologues) <- 
  c("Orthogroup", "Flybase_gene_ID", "Bter_gene_ID")

# FUNCTION DEFINITIONS ----

ReportPresenceOfAntioxidantGenes <- function (x, y) {
  # Reports whether a list of DEGs contains orthologues of the Dmel enzymatic 
  # antioxidant network genes.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting which treatment the B. terrestris DEG list is from
  #      ("C" or "R").
  #
  # Returns:
  #   Data.frame of results.
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Load DEGs list based on arguments.
  
  upregList <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                "up_regulated_LFC0.csv"))
  
  downregList <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                  "down_regulated_LFC0.csv"))
  
  # Add orthologues to the antioxidant gene list.
  
  antioxidantListOrthos <- merge(antioxidantGenes, dmelBterOrthologues)
  # 16 single-copy orthologues from the Dmel genes.
  
  # Loop over antioxidantListOrthos and add "1" to Study_expression column.
  
  for(i in (1:length(antioxidantListOrthos[, 1]))) {
    if (antioxidantListOrthos[i, "Bter_gene_ID"] %in% upregList[, "symbol"]) {
      antioxidantListOrthos[i, "Study_expression"] <- 1
    } else if (antioxidantListOrthos[i, "Bter_gene_ID"] %in% downregList[, "symbol"]) {
      antioxidantListOrthos[i, "Study_expression"] <- -1
    } else {
      antioxidantListOrthos[i, "Study_expression"] <- 0
    }
  }
  
  # Return antioxidantListOrthos.
  
  return(antioxidantListOrthos)
  
}

# EXECUTED STATEMENTS ----

# Bter Current Study ----

# Ovaries.

# Change working directory to 02_outputs.

setwd("../../../01_ovaries/31_DESeq2_DEG_lists/")

# Run function.

ovariesControlResults <- ReportPresenceOfAntioxidantGenes ("ovaries", "C")

ovariesRemovalResults <- ReportPresenceOfAntioxidantGenes ("ovaries", "R")

# Fat body.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

fatBodyControlResults <- ReportPresenceOfAntioxidantGenes ("fatbody", "C")

fatBodyRemovalResults <- ReportPresenceOfAntioxidantGenes ("fatbody", "R")

# Brain.

setwd("../../02_brain/31_DESeq2_DEG_lists/")

brainControlResults <- ReportPresenceOfAntioxidantGenes ("brain", "C")

brainRemovalResults <- ReportPresenceOfAntioxidantGenes ("brain", "R")

# Combine results.

combinedResults <- brainRemovalResults

colnames(combinedResults) <- 
  c("Flybase_gene_ID", "name", "Orthogroup", "Bter_gene_ID", "brain_removal")

# Add results.

combinedResults$brain_control <- brainControlResults[, 5]

combinedResults$fatbody_removal <- fatBodyRemovalResults[, 5]

combinedResults$fatbody_control <- fatBodyControlResults[, 5]

combinedResults$ovaries_removal <- ovariesRemovalResults[, 5]

combinedResults$ovaries_control <- ovariesControlResults[, 5]

# Plot and Save Heatmap ----

# Prepare data for heatmap.

# Heatmap showing genes which are differentially expressed in Bter.
# i.e. remove rows where all the DEG lists show 0 for expression.

noBterZeroResults <- combinedResults[rowSums(combinedResults[, 5:10] == 0, na.rm = TRUE) < 6, ]

# Plotting matrix of +1/0/-1 values.

plotMatrix <- as.matrix (noBterZeroResults[, c(5:10)])

# Add gene names to the matrix.

row.names(plotMatrix) <- noBterZeroResults[, "name"]

# Plot heatmap.

heatmapPlot <- PlotCustomHeatmap(46)

# Save heatmap.

setwd("../../00_all_tissues")

ggsave("38_NER0008751_obj1_fig_6b_antioxidant_heatmap.svg", plot = heatmapPlot, 
       width = 21, height = 29, units = "cm")
