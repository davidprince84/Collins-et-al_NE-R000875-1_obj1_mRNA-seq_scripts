#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Compare the differentially expressed genes (DEGs) between all samples  
# and samples with or without large numbers of virus-aligning reads.
#-------------------------------------------------------------------------------
# Inputs:
# Lists of DEGs from the current study.

# Outputs:
# .csv file summarising the results of the comparison in a table.
# .csv file recording all the DEGs that overlap between the two treatments.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# Data loaded by CompareVirusDEGs function.

# FUNCTION DEFINITIONS ----

CompareVirusDEGs <- function(x, y, z, a) {
  # Compares the overlap in DEG lists between all samples, no virus or with virus
  # samples for a given tissue, comparison of time points and direction of 
  # differential expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting which sets of DEGs to compare ("all_samples_vs_no_virus",
  #      ("all_samples_vs_with_virus", "no_virus_vs_with_virus")).
  #   z: string denoting the comparison between time points to be compared 
  #      ("R_TP1G vs. R_TP2G" or "R_TP2G vs. C_TP2G".)
  #   a: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in the earlier chronological time 
  #      point compared to the later chronological time point) or down-regulated
  #      genes (genes more highly expressed in the later chronological time point
  #      compared to the earlier chronological time point) ("up" or "down").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the comparison.
  #   2) Data.frame recording all the genes that overlap between the treatments. 
  # NOTE: The working directory needs to be set to 02_outputs/0x_tissue/33_DESeq2_additional_DEG_list 
  #       corresponding to the tissue used for argument x.
  
  # Set variables based on arguments.
  
  if (x == "brain") {
    supFileLetter <- "r"
  } else if (x == "fatbody") {
    supFileLetter <- "s"
  } else if (x == "ovaries") {
    supFileLetter <- "t"
  } else {
    stop ('Argument x must equal "brain", "fatbody" or "ovaries".')
  }
  
  
  if (y == "all_samples_vs_no_virus") {
    list1FileName <- paste0("00_NER0008751_obj1_sup_file_S1", supFileLetter, "_",
                            x, "_all_samples_additional_DEG_list.csv")
    list2FileName <- paste0("01_NER0008751_obj1_", x, "_no_virus_DEG_list.csv")
  } else if (y == "all_samples_vs_with_virus") {
    list1FileName <- paste0("00_NER0008751_obj1_sup_file_S1", supFileLetter, "_",
                            x, "_all_samples_additional_DEG_list.csv")
    list2FileName <- paste0("02_NER0008751_obj1_", x, "_with_virus_DEG_list.csv")
  } else if (y == "no_virus_vs_with_virus") {
    list1FileName <- paste0("01_NER0008751_obj1_", x, "_no_virus_DEG_list.csv")
    list2FileName <- paste0("02_NER0008751_obj1_", x, "_with_virus_DEG_list.csv")
  } else {
    stop('Argument y must equal "all_samples_vs_no_virus",
         "all_samples_vs_with_virus" or "no_virus_vs_with_virus".')
  }
  
  # Load DEG lists.
  
  list1DEGs <- read.csv(list1FileName)
  
  list2DEGs <- read.csv(list2FileName)
  
  # Index out DEGs for specific comparison and direction of differential
  # expression.
  
  compList1DEGs <- list1DEGs[list1DEGs[, "comparison"] == z, ]
  
  specificList1DEGs <- 
    compList1DEGs[compList1DEGs[, "expression_with_respect_to_chrono_age"] == paste0(a, "-regulated"), ]
  
  compList2DEGs <- list2DEGs[list2DEGs[, "comparison"] == z, ]
  
  specificList2DEGs <- 
    compList2DEGs[compList2DEGs[, "expression_with_respect_to_chrono_age"] == paste0(a, "-regulated"), ]
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Sets_of_DEGs" = y,
                          "Comparison" = z,
                          "Direction_of_expression" = paste0(a, "-regulated"),
                          "Number_of_list_1_DEGs" = length(specificList1DEGs$symbol),
                          "Number_of_list_2_DEGs" = length(specificList2DEGs$symbol),
                          "Number_overlapping_DEGs" = NA,
                          "Percentage_of_list_2_DEGs_overlapping" = NA) 
  
  # Calculate overlap between DEG lists.
  
  # Determine overlapping genes.
  
  overlappingGenes <- specificList2DEGs[specificList2DEGs$symbol %in% specificList1DEGs$symbol, ]
 
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_DEGs <- length(overlappingGenes$symbol)
  
  resultsDF$Percentage_of_list_2_DEGs_overlapping <- 
    round((resultsDF$Number_overlapping_DEGs/resultsDF$Number_of_list_2_DEGs)*100, digits = 1)
  
  # Prepare overlapping gene IDs for output.
  # Write results to data.frame, with extra blank space so that
  # all results can be combined later.
  
  genesToOutput <- overlappingGenes[, c("symbol", "name")]
  
  colnames(genesToOutput) <- c(paste(x, y, z, a, "regulated_Bter_gene_symbol", sep = "_"),
                               paste(x, y, z, a, "regulated_Bter_gene_name", sep = "_"))
  
  if (is.na(genesToOutput[1, 1])) {
    genesToOutput[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(genesToOutput[, 1]) + 1
  
  genesToOutput[fillStart:4000, ] <- ""
  
  # Make and return the list of outputs.
  
  listToReturn <- list("genes" = genesToOutput, 
                       "comparison" = resultsDF)
  
  return(listToReturn)
  
}

# EXECUTED STATEMENTS ----

# Brain DEG Overlap ----

# Change working directory.
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

setwd("../02_outputs/02_brain/33_DESeq2_additional_DEG_list/")

# Compare DEG lists.

brainAllvsNoVirusRTP1vsRTP2UpRegOverlap <- 
  CompareVirusDEGs("brain", "all_samples_vs_no_virus",
                    "R_TP1G vs. R_TP2G", "up")

brainAllvsNoVirusRTP1vsRTP2DownRegOverlap <- 
  CompareVirusDEGs("brain", "all_samples_vs_no_virus",
                    "R_TP1G vs. R_TP2G", "down")

brainAllvsNoVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("brain", "all_samples_vs_no_virus",
                   "R_TP2G vs. C_TP2G", "up")

brainAllvsNoVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("brain", "all_samples_vs_no_virus",
                   "R_TP2G vs. C_TP2G", "down")

brainAllvsWithVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("brain", "all_samples_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "up")

brainAllvsWithVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("brain", "all_samples_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "down")

brainNoVirusvsWithVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("brain", "no_virus_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "up")

brainNoVirusvsWithVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("brain", "no_virus_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "down")

# Fat Body DEG Overlap ----

# Change directory.

setwd("../../03_fatbody/33_DESeq2_additional_DEG_list/")

# Compare DEG lists.

fatBodyAllvsNoVirusRTP1vsRTP2UpRegOverlap <- 
  CompareVirusDEGs("fatbody", "all_samples_vs_no_virus",
                   "R_TP1G vs. R_TP2G", "up")

fatBodyAllvsNoVirusRTP1vsRTP2DownRegOverlap <- 
  CompareVirusDEGs("fatbody", "all_samples_vs_no_virus",
                   "R_TP1G vs. R_TP2G", "down")

fatBodyAllvsNoVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("fatbody", "all_samples_vs_no_virus",
                   "R_TP2G vs. C_TP2G", "up")

fatBodyAllvsNoVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("fatbody", "all_samples_vs_no_virus",
                   "R_TP2G vs. C_TP2G", "down")

# Ovaries DEG Overlap ----

# Change directory.

setwd("../../01_ovaries/33_DESeq2_additional_DEG_list/")

# Compare DEG lists.

ovariesAllvsNoVirusRTP1vsRTP2UpRegOverlap <- 
  CompareVirusDEGs("ovaries", "all_samples_vs_no_virus", "R_TP1G vs. R_TP2G", "up")

ovariesAllvsNoVirusRTP1vsRTP2DownRegOverlap <- 
  CompareVirusDEGs("ovaries", "all_samples_vs_no_virus", "R_TP1G vs. R_TP2G", "down")

ovariesAllvsNoVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("ovaries", "all_samples_vs_no_virus", "R_TP2G vs. C_TP2G", "up")

ovariesAllvsNoVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("ovaries", "all_samples_vs_no_virus", "R_TP2G vs. C_TP2G", "down")

ovariesAllvsWithVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("ovaries", "all_samples_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "up")

ovariesAllvsWithVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("ovaries", "all_samples_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "down")

ovariesNoVirusvsWithVirusRTP2vsCTP2UpRegOverlap <- 
  CompareVirusDEGs("ovaries", "no_virus_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "up")

ovariesNoVirusvsWithVirusRTP2vsCTP2DownRegOverlap <- 
  CompareVirusDEGs("ovaries", "no_virus_vs_with_virus",
                   "R_TP2G vs. C_TP2G", "down")

# Combine Results into Single Data.Frames ----

# Combine stats results.

combinedStats <- rbind(brainAllvsNoVirusRTP1vsRTP2UpRegOverlap$comparison,
                       brainAllvsNoVirusRTP1vsRTP2DownRegOverlap$comparison,
                       brainAllvsNoVirusRTP2vsCTP2UpRegOverlap$comparison,
                       brainAllvsNoVirusRTP2vsCTP2DownRegOverlap$comparison,
                       brainAllvsWithVirusRTP2vsCTP2UpRegOverlap$comparison,
                       brainAllvsWithVirusRTP2vsCTP2DownRegOverlap$comparison,
                       brainNoVirusvsWithVirusRTP2vsCTP2UpRegOverlap$comparison,
                       brainNoVirusvsWithVirusRTP2vsCTP2DownRegOverlap$comparison,
                       fatBodyAllvsNoVirusRTP1vsRTP2UpRegOverlap$comparison,
                       fatBodyAllvsNoVirusRTP1vsRTP2DownRegOverlap$comparison,
                       fatBodyAllvsNoVirusRTP2vsCTP2UpRegOverlap$comparison,
                       fatBodyAllvsNoVirusRTP2vsCTP2DownRegOverlap$comparison,
                       ovariesAllvsNoVirusRTP1vsRTP2UpRegOverlap$comparison,
                       ovariesAllvsNoVirusRTP1vsRTP2DownRegOverlap$comparison,
                       ovariesAllvsNoVirusRTP2vsCTP2UpRegOverlap$comparison,
                       ovariesAllvsNoVirusRTP2vsCTP2DownRegOverlap$comparison,
                       ovariesAllvsWithVirusRTP2vsCTP2UpRegOverlap$comparison,
                       ovariesAllvsWithVirusRTP2vsCTP2DownRegOverlap$comparison,
                       ovariesNoVirusvsWithVirusRTP2vsCTP2UpRegOverlap$comparison,
                       ovariesNoVirusvsWithVirusRTP2vsCTP2DownRegOverlap$comparison)

# Combine overlapping genes results.

combinedGenes <- cbind(brainAllvsNoVirusRTP1vsRTP2UpRegOverlap$genes,
                       brainAllvsNoVirusRTP1vsRTP2DownRegOverlap$genes,
                       brainAllvsNoVirusRTP2vsCTP2UpRegOverlap$genes,
                       brainAllvsNoVirusRTP2vsCTP2DownRegOverlap$genes,
                       brainAllvsWithVirusRTP2vsCTP2UpRegOverlap$genes,
                       brainAllvsWithVirusRTP2vsCTP2DownRegOverlap$genes,
                       brainNoVirusvsWithVirusRTP2vsCTP2UpRegOverlap$genes,
                       brainNoVirusvsWithVirusRTP2vsCTP2DownRegOverlap$genes,
                       fatBodyAllvsNoVirusRTP1vsRTP2UpRegOverlap$genes,
                       fatBodyAllvsNoVirusRTP1vsRTP2DownRegOverlap$genes,
                       fatBodyAllvsNoVirusRTP2vsCTP2UpRegOverlap$genes,
                       fatBodyAllvsNoVirusRTP2vsCTP2DownRegOverlap$genes,
                       ovariesAllvsNoVirusRTP1vsRTP2UpRegOverlap$genes,
                       ovariesAllvsNoVirusRTP1vsRTP2DownRegOverlap$genes,
                       ovariesAllvsNoVirusRTP2vsCTP2UpRegOverlap$genes,
                       ovariesAllvsNoVirusRTP2vsCTP2DownRegOverlap$genes,
                       ovariesAllvsWithVirusRTP2vsCTP2UpRegOverlap$genes,
                       ovariesAllvsWithVirusRTP2vsCTP2DownRegOverlap$genes,
                       ovariesNoVirusvsWithVirusRTP2vsCTP2UpRegOverlap$genes,
                       ovariesNoVirusvsWithVirusRTP2vsCTP2DownRegOverlap$genes)

# Save Results ----

setwd("../../00_all_tissues/")

write.csv(combinedStats, "25_NER0008751_obj1_all_samples_vs_with_without_virus_comparison_stats.csv", 
          row.names = FALSE)

write.csv(combinedGenes, "26_NER0008751_obj1_all_samples_vs_with_without_virus_overlapping_genes.csv", 
          row.names = FALSE)
