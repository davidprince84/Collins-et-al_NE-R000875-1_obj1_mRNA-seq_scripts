#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Comparison of Bombus terrestris (Bter) differentially expressed genes 
# (DEGs) with Drosophila melanogaster (Dmel) differentially expressed genes from
# Chen et al. (2014) and Pacifico et al. (2018).
#-------------------------------------------------------------------------------
# Inputs:
# OrthoFinder results. Lists of Dmel DEGs from Chen et al. (2014) and 
# Pacifico et al. (2018). Lists of Bter DEGs from the current study.

# Outputs:
# .csv file summarising the results of the statistical tests in a table.
# .csv file recording all the DEGs that overlap between Dmel and Bter.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.
# org.Dm.eg.db loaded via the script in the "Load data" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("40_compare_dmel_genage_prep.R")

# Format Objects as Required ----

colnames(dmelBterOrthologues) <- c("Orthogroup",
                                   "Dmel_gene_ID",
                                   "symbol")

# FUNCTION DEFINITIONS ----

CompareBterAndDmelDEGs <- function (x, y, z) {
  # Compares the overlap in DEG lists between the Bter data from the current 
  # study and Dmel data from previous studies for a given tissue, and direction 
  # of differential expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain" or
  #      "fatbody").
  #   y: string denoting which treatment the B. terrestris DEG list is from
  #      ("C" or "R").
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in time point 2 compared to time 
  #      point one) or down-regulated genes (genes more highly expressed in time
  #      point 1 compared to time point 2) ("up" or "down").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between Dmel and Bter. 
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Set variables based on arguments.
  
  if (x == "fatbody" && z == "up") {
    dmelList <- chenUpDEGsFB
  } else if (x == "fatbody" && z == "down") {
    dmelList <- chenDownDEGsFB 
  } else if (x == "brain" && z == "up") {
    dmelList <- pacificoUpDEGsFB
  } else if (x == "brain" && z == "down") {
    dmelList <- pacificoDownDEGsFB
  } else {
    stop ("Argument x or z is incorrect.")
  }
  
  # Load Bter lists.
  
  bterToCompare <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                   z , "_regulated_LFC0.csv"))
  
  backgroundList <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Reduce number of columns and rename.
  
  backgroundList <- as.data.frame(backgroundList[, 1])
  
  colnames(backgroundList) <- "symbol"
  
  # Add Dmel orthologues.
  
  bterWithOrthos <- merge(bterToCompare, dmelBterOrthologues)
  
  backgroundWithOrthos <- merge(backgroundList, dmelBterOrthologues)
  
  # Dmel list to compare.
  
  dmelListOrthos <- merge(dmelList, dmelBterOrthologues)
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Treatment" = y,
                          "Direction_of_expression" = paste0(z, "-regulated"),
                          "Number_Bter_DEGs" = length(bterWithOrthos$symbol),
                          "Number_Dmel_DEGs" = length(dmelListOrthos$symbol),
                          "Number_overlapping_DEGs" = NA,
                          "Number_DEGs_only_in_Bter" = NA,
                          "Number_DEGs_only_in_Dmel" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = (0.05/4),
                          "Percentage_of_Bter_DEGs_overlapping" = NA) 
  
  # Calculate overlap between Bter and Dmel.
  
  # Determine overlapping genes.
  
  overlappingGenes <- bterWithOrthos[bterWithOrthos$symbol %in% dmelListOrthos$symbol, ]
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$symbol)
  
  resultsDF$Percentage_of_Bter_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_Bter_DEGs)*100, digits = 1)
  
  # Prepare overlapping gene IDs for output.
  # Write results to data.frame, with extra blank space so that
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes[, c("symbol", "name", "Dmel_gene_ID")]
  
  colnames(genesToOutput) <- c(paste(x, y, z, "regulated_Bter_gene_ID", sep = "_"),
                               paste(x, y, z, "regulated_Bter_gene_name", sep = "_"),
                               paste(x, y, z, "regulated_Flybase_gene_ID", sep = "_"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:200, ] <- ""
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_DEGs_only_in_Bter <- 
    resultsDF$Number_Bter_DEGs - resultsDF$Number_overlapping_DEGs
  
  resultsDF$Number_DEGs_only_in_Dmel <-
    resultsDF$Number_Dmel_DEGs - resultsDF$Number_overlapping_DEGs
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundWithOrthos$symbol) - resultsDF$Number_overlapping_DEGs - 
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
  
  # Make and return the list of outputs.
  
  listToReturn <- list("genes" = genesToOutput, 
                       "stats" = resultsDF)
  
  return(listToReturn)
  
}

# EXECUTED STATEMENTS ----

# Brain DEG Overlap ----

# Change working directory.

setwd("../02_outputs/02_brain/31_DESeq2_DEG_lists/")

# Compare DEG lists.

brainControlUpOverlap <- CompareBterAndDmelDEGs("brain", "C", "up")

brainControlDownOverlap <- CompareBterAndDmelDEGs("brain", "C", "down")

brainRemovalUpOverlap <- CompareBterAndDmelDEGs("brain", "R", "up")

brainRemovalDownOverlap <- CompareBterAndDmelDEGs("brain", "R", "down")

# Fat Body DEG Overlap ----

# Change directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare DEG lists.

fatBodyControlUpOverlap <- CompareBterAndDmelDEGs("fatbody", "C", "up")

fatBodyControlDownOverlap <- CompareBterAndDmelDEGs("fatbody", "C", "down")

fatBodyRemovalUpOverlap <- CompareBterAndDmelDEGs("fatbody", "R", "up")

fatBodyRemovalDownOverlap <- CompareBterAndDmelDEGs("fatbody", "R", "down")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(brainControlUpOverlap$stats,
                       brainControlDownOverlap$stats,
                       brainRemovalUpOverlap$stats,
                       brainRemovalDownOverlap$stats,
                       fatBodyControlUpOverlap$stats,
                       fatBodyControlDownOverlap$stats,
                       fatBodyRemovalUpOverlap$stats,
                       fatBodyRemovalDownOverlap$stats)

# Combine genes results.

combinedGenes <- cbind(brainControlUpOverlap$genes,
                       brainControlDownOverlap$genes,
                       brainRemovalUpOverlap$genes,
                       brainRemovalDownOverlap$genes,
                       fatBodyControlUpOverlap$genes,
                       fatBodyControlDownOverlap$genes,
                       fatBodyRemovalUpOverlap$genes,
                       fatBodyRemovalDownOverlap$genes)

# Save Results ----

setwd("../../00_all_tissues/")

write.csv(combinedStats, "39_NER0008751_obj1_table_S10_Dmel_stats_results.csv", 
          row.names = FALSE)

write.csv(combinedGenes, "40_NER0008751_obj1_table_S11_overlapping_genes.csv", 
          row.names = FALSE)
