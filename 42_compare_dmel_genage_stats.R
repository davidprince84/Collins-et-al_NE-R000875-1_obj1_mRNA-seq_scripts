#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Statistical comparison of overlap between Bombus terrestris (Bter) 
# differentially expressed genes (DEGs) with Drosophila melanogaster (Dmel) 
# GenAge data set.
#-------------------------------------------------------------------------------
# Inputs:
# GenAge Dmel gene list with Bter orthologs (from script #40). Lists of 
# Bter DEGs.

# Outputs:
# .csv file summarising the results of the statistical tests in a table.
# .csv file recording all the DEGs that overlap between the Dmel and Bter.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# org.Dm.eg.db loaded via the script in the "Load data" section.

# LOAD DATA ----

# Load GenAge Data ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("40_compare_dmel_genage_prep.R")

# Remove unnecessary data.frames.

rm(chenDownDEGsFB, chenUpDEGsFB, pacificoUpDEGsFB, pacificoDownDEGsFB)

# FUNCTION DEFINITIONS ----

CompareBterDEGsAndDmelGenAge <- function (x, y) {
  # Compares the overlap in single-copy orthologous genes between Bter DEG lists
  # from the current study and Dmel genes in the GenAge database.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting which treatment the B. terrestris DEG list is from
  #      ("C" or "R").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between Dmel and Bter. 
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Load Bter lists.
  
  bterUpDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                "up_regulated_LFC0.csv"))
  
  bterDownDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                  "down_regulated_LFC0.csv"))
  
  backgroundList <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Reduce number of columns and rename column in background list.
  
  backgroundList <- as.data.frame(backgroundList[, 1])
  
  colnames(backgroundList) <- "Bter_gene_ID"
  
  # Combine Bter up- and down-regulated genes, simplify and rename column.
  
  bterDEGs <- rbind(bterUpDEGs, bterDownDEGs)
  
  bterDEGs <- as.data.frame(bterDEGs[, c("symbol", "name")])
  
  colnames(bterDEGs) <- c("Bter_gene_ID", "Bter_gene_name")
  
  # Add Dmel orthologues.
  
  bterWithOrthos <- merge(bterDEGs, dmelBterOrthologues)
  
  backgroundWithOrthos <- merge(backgroundList, dmelBterOrthologues)
  
  # Dmel list to compare.
  
  dmelListOrthos <- genAgeOrtho
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Treatment" = y,
                          "Dmel_list" = "GenAge",
                          "Number_Bter_DEGs" = length(bterWithOrthos$Bter_gene_ID),
                          "Number_Dmel_DEGs" = length(dmelListOrthos$Flybase_gene_ID),
                          "Number_overlapping_DEGs" = NA,
                          "Number_DEGs_only_in_Bter" = NA,
                          "Number_DEGs_only_in_Dmel" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = (0.05/6),
                          "Percentage_of_Bter_DEGs_overlapping" = NA) 
  
  # Calculate overlap between Bter and Dmel.
  
  # Determine overlapping genes.
  
  overlappingGenes <- bterWithOrthos[bterWithOrthos$Bter_gene_ID %in% dmelListOrthos$Bter_gene_ID, ]
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$Bter_gene_ID)
  
  resultsDF$Percentage_of_Bter_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_Bter_DEGs)*100, digits = 1)
  
  # Write results to data.frame, with extra blank space so that 
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes[, c(1, 2, 4)]
  
  colnames(genesToOutput) <- c(paste(x, y, "Bter_gene_ID", sep = "_"),
                               paste(x, y, "Bter_gene_name", sep = "_"),
                               paste(x, y, "Flybase_gene_ID", sep = "_"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:length(dmelListOrthos[, 1]), ] <- ""
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_DEGs_only_in_Bter <- 
    resultsDF$Number_Bter_DEGs - resultsDF$Number_overlapping_DEGs
  
  resultsDF$Number_DEGs_only_in_Dmel <-
    resultsDF$Number_Dmel_DEGs - resultsDF$Number_overlapping_DEGs
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundWithOrthos$Bter_gene_ID) - resultsDF$Number_overlapping_DEGs - 
    resultsDF$Number_DEGs_only_in_Bter - resultsDF$Number_DEGs_only_in_Dmel
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both DEG lists, only 
  # control DEG list, only removal DEG list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_DEGs, 
                               resultsDF$Number_DEGs_only_in_Bter,
                               resultsDF$Number_DEGs_only_in_Dmel, 
                               resultsDF$Number_genes_in_neither))
  
  # Change the dimensions to 2 rows and 2 columns.
  
  dim(contingencyTable) <- c(2,2)
  
  # Conduct two-tailed Fisher's Exact Test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance.
  
  fisherResults <- fisher.test(contingencyTable)
  
  # Add results of statistical tests to resultsDF.
  
  resultsDF$p_value <- fisherResults$p.value
  
  resultsDF$Odds_ratio <- fisherResults$estimate[[1]]
  
  # Return list of results.
  
  listToReturn <- list("genes" = genesToOutput, 
                       "stats" = resultsDF)
  
  return(listToReturn)
  
}

# EXECUTED STATEMENTS ----

# Brain Comparison ----

# Set directory.

setwd("../02_outputs/02_brain/31_DESeq2_DEG_lists/")

# Compare gene lists.

brainControlResults <- CompareBterDEGsAndDmelGenAge("brain", "C")

brainRemovalResults <- CompareBterDEGsAndDmelGenAge("brain", "R")

# Fat Body Comparison ----

# Set directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare gene lists.

fatBodyControlResults <- CompareBterDEGsAndDmelGenAge("fatbody", "C")

fatBodyRemovalResults <- CompareBterDEGsAndDmelGenAge("fatbody", "R")

# Ovaries Comparison ----

# Set directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists/")

# Compare gene lists.

ovariesControlResults <- CompareBterDEGsAndDmelGenAge("ovaries", "C")

ovariesRemovalResults <- CompareBterDEGsAndDmelGenAge("ovaries", "R")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(brainControlResults$stats, brainRemovalResults$stats,
                       fatBodyControlResults$stats, fatBodyRemovalResults$stats,
                       ovariesControlResults$stats, ovariesRemovalResults$stats)

# Combine genes results.

combinedGenes <- cbind(brainControlResults$genes, brainRemovalResults$genes,
                       fatBodyControlResults$genes, fatBodyRemovalResults$genes,
                       ovariesControlResults$genes, ovariesRemovalResults$genes)

# Save Results ----

setwd("../../00_all_tissues/")

write.csv(combinedStats, "31_NER0008751_obj1_table_S12_GenAge_stats.csv",
          row.names = FALSE)

write.csv(combinedGenes, "32_NER0008751_obj1_table_S13_GenAge_overlapping_genes.csv",
          row.names = FALSE)
