#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Compare the differentially expressed genes (DEGs) between chronological
# ages within tissues.
#-------------------------------------------------------------------------------
# Inputs:
# Lists of DEGs from the current study.

# Outputs:
# .csv file summarising the results of the statistical tests in a table.
# .csv file recording all the DEGs that overlap between the two treatments.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# Data loaded by CompareChronoAgeDEGs function.

# FUNCTION DEFINITIONS ----

CompareChronoAgeDEGs <- function (x, y, z) {
  # Compares the overlap in DEG lists between a within treatment comparison list
  # and a between treatment comparison list with similar chronological age 
  # differences for a given tissue, and direction of differential expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in the later time-point compared to 
  #      the earlier time-point) or down-regulated genes (genes more highly 
  #      expressed in the later time-point compared to the earlier time-point)
  #      ("up" or "down").
  #   z: string denoting whether the DEG lists to compare are for the early 
  #      section of chronological age (R:TP1-C:TP1 and R:TP1-R:TP2) or the late 
  #      section of chronological age (R:TP2-C:TP2 and C:TP1-C:TP2) ("early" or
  #      "late").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between the treatments. 
  # NOTE: The working directory needs to be set output foldee for the specified 
  # tissue (e.g. 02_outputs/01_ovaries for ovaries).
  
  # Set variables based on arguments x, y and z.
  
  if (x == "brain") {
    supFileLetter <- "r"
  } else if (x == "fatbody") {
    supFileLetter <- "s"
  } else if (x == "ovaries") {
    supFileLetter <- "t"
  } else {
    stop ('Argument x must equal "brain", "fatbody" or "ovaries".')
  }
  
  if (y == "up") {
    directionIndex <- "up-regulated"
  } else if (y == "down") {
    directionIndex <- "down-regulated"
  } else {
    stop ('Argument y must equal "up" or "down".')
  }
  
  if (z == "early") {
    treatment <- "R"
    compIndex <- "R_TP1G vs. C_TP1G"
    withinName <- "R_TP1G vs. R_TP2G"
  } else if (z == "late") {
    treatment <- "C"
    compIndex <- "R_TP2G vs. C_TP2G"
    withinName <- "C_TP1G vs. C_TP2G"
  } else {
    stop ('Argument z must equal "early" or "late".')
  }
  
  # Load data based on arguments.
  
  # Individual within treatment DEGs.
  
  withinTreatmentDEGs <- read.csv(paste0("31_DESeq2_DEG_lists/", x, "_results_", 
                                         treatment, "_treatment_", y,
                                         "_regulated_LFC0.csv"))
  
  # Index between treatment DEGs from a combined file for the tissue.
  
  allBetweenTreatmentDEGs <- read.csv(paste0("33_DESeq2_additional_DEG_list/",
                                             "00_NER0008751_obj1_sup_file_S1", 
                                             supFileLetter, "_",
                                             x, "_all_samples_additional_DEG_list.csv"))
  
  compBetweenTreatmentDEGs <- 
    allBetweenTreatmentDEGs[allBetweenTreatmentDEGs[, "comparison"] == compIndex, ]
  
  betweenTreatmentDEGs <- 
    compBetweenTreatmentDEGs[compBetweenTreatmentDEGs[, "expression_with_respect_to_chrono_age"] == directionIndex, ]
  
  # Background list of expressed genes for statistical test.
  
  backgroundList <- read.csv(paste0("31_DESeq2_DEG_lists/",
                                    x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Initialize results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Within_treatment_DEGs" = withinName,
                          "Between_treatment_DEGs" = compIndex,
                          "Direction_of_expression" = directionIndex,
                          "Number_within_DEGs" = length(withinTreatmentDEGs$symbol),
                          "Number_between_DEGs" = length(betweenTreatmentDEGs$symbol),
                          "Number_overlapping_DEGs" = NA,
                          "Number_DEGs_only_in_within" = NA,
                          "Number_DEGs_only_in_between" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = (0.05/4),
                          "Percentage_of_within_DEGs_overlapping" = NA) 
  
  # Calculate overlap between control and removal treatments.
  
  # Determine overlapping genes.
  
  overlappingGenes <- 
    withinTreatmentDEGs[withinTreatmentDEGs$symbol %in% betweenTreatmentDEGs$symbol, ]
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$symbol)
  
  resultsDF$Percentage_of_within_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_within_DEGs)*100, digits = 1)
  
  # Prepare overlapping gene IDs for output.
  # Write results to data.frame, with extra blank space so that
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes[, c("symbol", "name")]
  
  colnames(genesToOutput) <- c(paste(x, z, y, "regulated_Bter_gene_symbol", sep = "_"),
                               paste(x, z, y, "regulated_Bter_gene_name", sep = "_"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:4000, ] <- ""
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_DEGs_only_in_within <- 
    resultsDF$Number_within_DEGs - resultsDF$Number_overlapping_DEGs
  
  resultsDF$Number_DEGs_only_in_between <-
    resultsDF$Number_between_DEGs - resultsDF$Number_overlapping_DEGs
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundList$gene_symbol) - resultsDF$Number_overlapping_DEGs - 
    resultsDF$Number_DEGs_only_in_within - resultsDF$Number_DEGs_only_in_between
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both DEG lists, only 
  # within DEG list, only between DEG list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_DEGs, 
                               resultsDF$Number_DEGs_only_in_within,
                               resultsDF$Number_DEGs_only_in_between, 
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
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

setwd("../02_outputs/02_brain/")

# Compare DEG lists.

brainEarlyUpRegOverlap <- CompareChronoAgeDEGs("brain", "up", "early")

brainEarlyDownRegOverlap <- CompareChronoAgeDEGs("brain", "down", "early")

brainLateUpRegOverlap <- CompareChronoAgeDEGs("brain", "up", "late")

brainLateDownRegOverlap <- CompareChronoAgeDEGs("brain", "down", "late")

# Fat Body DEG Overlap ----

# Change directory.

setwd("../03_fatbody/")

# Compare DEG lists.

fatBodyEarlyUpRegOverlap <- CompareChronoAgeDEGs("fatbody", "up", "early")

fatBodyEarlyDownRegOverlap <- CompareChronoAgeDEGs("fatbody", "down", "early")

fatBodyLateUpRegOverlap <- CompareChronoAgeDEGs("fatbody", "up", "late")

fatBodyLateDownRegOverlap <- CompareChronoAgeDEGs("fatbody", "down", "late")

# Ovaries DEG Overlap ----

# Change directory.

setwd("../01_ovaries/")

# Compare DEG lists.

ovariesEarlyUpRegOverlap <- CompareChronoAgeDEGs("ovaries", "up", "early")

ovariesEarlyDownRegOverlap <- CompareChronoAgeDEGs("ovaries", "down", "early")

ovariesLateUpRegOverlap <- CompareChronoAgeDEGs("ovaries", "up", "late")

ovariesLateDownRegOverlap <- CompareChronoAgeDEGs("ovaries", "down", "late")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(brainEarlyUpRegOverlap$stats,
                       brainEarlyDownRegOverlap$stats,
                       brainLateUpRegOverlap$stats,
                       brainLateDownRegOverlap$stats,
                       fatBodyEarlyUpRegOverlap$stats,
                       fatBodyEarlyDownRegOverlap$stats,
                       fatBodyLateUpRegOverlap$stats,
                       fatBodyLateDownRegOverlap$stats,
                       ovariesEarlyUpRegOverlap$stats,
                       ovariesEarlyDownRegOverlap$stats,
                       ovariesLateUpRegOverlap$stats,
                       ovariesLateDownRegOverlap$stats)

# Combine overlapping genes results.

combinedGenes <- cbind(brainEarlyUpRegOverlap$genes,
                       brainEarlyDownRegOverlap$genes,
                       brainLateUpRegOverlap$genes,
                       brainLateDownRegOverlap$genes,
                       fatBodyEarlyUpRegOverlap$genes,
                       fatBodyEarlyDownRegOverlap$genes,
                       fatBodyLateUpRegOverlap$genes,
                       fatBodyLateDownRegOverlap$genes,
                       ovariesEarlyUpRegOverlap$genes,
                       ovariesEarlyDownRegOverlap$genes,
                       ovariesLateUpRegOverlap$genes,
                       ovariesLateDownRegOverlap$genes)

# Save Results ----

setwd("../00_all_tissues/")

write.csv(combinedStats, "28_NER0008751_obj1_table_S21_age_comparison_stats.csv", 
          row.names = FALSE)

write.csv(combinedGenes, "29_NER0008751_obj1_table_S22_age_overlapping_genes.csv", 
          row.names = FALSE)
