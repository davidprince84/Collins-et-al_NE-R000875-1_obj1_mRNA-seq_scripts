#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce lists of differentially expressed genes (DEGs) for each tissue
# comparing within treatments.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file, genesymbols2genenames.txt file, Kallisto 
# abundances, virus_samples.csv

# Outputs:
# .csv files containing lists of DEGs for a given tissue and treatment.
# .csv file containing all genes expressed in the analysis, to be used as a
# background list in downstream comparisons.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# DESeq2 and tximport loaded via the scripts in the 
# "Custom functions" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Transcript to Gene File ----

setwd("../00_data/10_transcripts2genes")

t2g <- read.table(file = "Bter_v1_transcripts2genes.txt",
                  header = FALSE,
                  col.names = c("TXNAME",
                                "GENEID"))

# Gene Symbol to Gene Name File ----

setwd("../05_feature_table/")

symbolToName <- read.csv("genesymbols2genenames.txt")

# Virus Samples ----

setwd("../11_virus_samples")

virus_samples <- read.csv(file = "virus_samples.csv",
                              col.names = c("sample", "virus"))

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_BuildDDSForDEAnalysis.R")

# FUNCTION DEFINITIONS ----

MakeGeneList <- function (x, LFC = 0) {
  # Produces and saves a list of all expressed genes for a tissue with 
  # differential expression results added.
  #
  # Args:
  #   x: string denoting which tissue to generate a DEG list for ("brain",
  #      "fatbody" or "ovaries").
  #   LFC: numeric denoting the LFC to use when filter DESeq2 results
  #      (default of 0 (i.e. all significant genes)).
  #
  # Returns:
  #   A .csv file of Bombus terrestris genes with associated normalized counts 
  #   for each sample and differential expression results.
  
  # Note: false discovery rate (FDR) of 0.05 used, but this can be changed by
  # changing the alpha level in the results() call below.
  
  # Set variables based on arguments.
  
  if (x == "brain") {
    dds <- ddsBrain
    tableNumber <- "S4"
  } else if (x == "fatbody") {
    dds <- ddsFatBody
    tableNumber <- "S5"
  } else if (x == "ovaries") {
    dds <- ddsOvaries
    tableNumber <- "S6"
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  # Generate data.frame of normalized counts for all genes and add gene names
  # as a column (for merger with differential expression results).
  
  normalizedCounts <- as.data.frame(counts(dds, normalized = TRUE))
  
  normalizedCounts$gene_symbol <- row.names(normalizedCounts)
  
  # Extract differential expression results for C (control) and add gene names
  # as a column (for merger with normalized count data).
  
  allCExtractedResults <- 
    as.data.frame(results(dds, contrast = c("condition", "C_TP2G", "C_TP1G"),
                          alpha = 0.05, lfcThreshold = LFC))
  
  allCExtractedResults$gene_symbol <- row.names(allCExtractedResults)
  
  # Add columns explicitly stating where a gene is differentially expressed and
  # whether it is up-regulated (with age, i.e. more expressed in old bees)
  # or down-regulated (with age, i.e. more expressed in younger bees).
  
  allCExtractedResults[, "DEG"] <- NA
  allCExtractedResults[, "expression_with_age"] <- NA
  
  for (i in (1:length(allCExtractedResults[, 1]))) {
    if (is.na(allCExtractedResults[i, "padj"])) {
      allCExtractedResults[i, "DEG"] <- "NO"
    } else if (allCExtractedResults[i, "padj"] < 0.05) {
      allCExtractedResults[i, "DEG"] <- "YES"
    } else {
      allCExtractedResults[i, "DEG"] <- "NO"
    }
    if (allCExtractedResults[i, "DEG"] == "YES" && allCExtractedResults[i, "log2FoldChange"] > 0) {
      allCExtractedResults[i, "expression_with_age"] <- "up-regulated"
    } else if (allCExtractedResults[i, "DEG"] == "YES" && allCExtractedResults[i, "log2FoldChange"] < 0) {
      allCExtractedResults[i, "expression_with_age"] <- "down-regulated"
    } else {
      allCExtractedResults[i, "expression_with_age"] <- "NA"
    }
  }
  
  # Rename columns.
  
  colnames(allCExtractedResults) <- 
    c("baseMean_C", "log2FoldChange_C", "lfcSE_C", "stat_C", "pvalue_C",
      "padj_C", "gene_symbol", "DEG_in_C", "expression_with_age_in_C")
  
  # Extract differential expression results for R (removal) and add gene names
  # as a column (for merger with normalized count data).
  
  allRExtractedResults <- 
    as.data.frame(results(dds, contrast = c("condition", "R_TP2G", "R_TP1G"),
                          alpha = 0.05, lfcThreshold = LFC))
  
  allRExtractedResults$gene_symbol <- row.names(allRExtractedResults)
  
  # Add columns explicitly stating where a gene is differentially expressed and
  # whether it is up-regulated (with age, i.e. more expressed in old bees)
  # or down-regulated (with age, i.e. more expressed in younger bees).
  
  allRExtractedResults[, "DEG"] <- NA
  allRExtractedResults[, "expression_with_age"] <- NA
  
  for (i in (1:length(allRExtractedResults[, 1]))) {
    if (is.na(allRExtractedResults[i, "padj"])) {
      allRExtractedResults[i, "DEG"] <- "NO"
    } else if (allRExtractedResults[i, "padj"] < 0.05) {
      allRExtractedResults[i, "DEG"] <- "YES"
    } else {
      allRExtractedResults[i, "DEG"] <- "NO"
    }
    if (allRExtractedResults[i, "DEG"] == "YES" && allRExtractedResults[i, "log2FoldChange"] > 0) {
      allRExtractedResults[i, "expression_with_age"] <- "up-regulated"
    } else if (allRExtractedResults[i, "DEG"] == "YES" && allRExtractedResults[i, "log2FoldChange"] < 0) {
      allRExtractedResults[i, "expression_with_age"] <- "down-regulated"
    } else {
      allRExtractedResults[i, "expression_with_age"] <- "NA"
    }
  }
  
  # Rename columns.
  
  colnames(allRExtractedResults) <- 
    c("baseMean_R", "log2FoldChange_R", "lfcSE_R", "stat_R", "pvalue_R",
      "padj_R", "gene_symbol", "DEG_in_R", "expression_with_age_in_R")
  
  # Merge normalizedCounts, allCExtractedResults and allRExtractedResults to 
  # attach count data to differential expression results.
  
  countsAndResults <- merge(normalizedCounts, allCExtractedResults, by = "gene_symbol")
  
  countsAndResults <- merge(countsAndResults, allRExtractedResults, by = "gene_symbol")
  
  # Merge countsAndResults with symbolToName to add gene names.
  
  colnames(symbolToName) <- c("gene_symbol", "gene_name")
  
  namedResults <- merge(countsAndResults, symbolToName, by.x = "gene_symbol",
                        all.x = TRUE)
  
  # Replace NAs in name with "uncharacterized $gene_symbol".
  
  for (row in 1:length(namedResults$gene_symbol)) {
    if (is.na(namedResults[row, "gene_name"])) {
        namedResults[row, "gene_name"] <- paste0("uncharacterized ", namedResults[row, "gene_symbol"])
    }
  }
  
  # Write namedResults to file for manuscript supplementary data and 
  # submission to Gene Expression Omnibus (GEO).
  
  write.csv(namedResults, 
            paste0(x, "_all_expressed_genes_and_DE_results.csv"),
            row.names = FALSE)
  
  write.csv(namedResults, 
            paste0("00_NER0008751_obj1_table_", tableNumber, ".csv"),
            row.names = FALSE)
}

MakeDEGLists <- function (x) {
  # Make DEG lists for a tissue by dividing up the list of all genes into
  # lists of DEGs.
  #
  # Args:
  #   x: string denoting which tissue to generate a DEG list for ("brain",
  #      "fatbody" or "ovaries").
  #
  # Returns:
  #   Four .csv files of Bombus terrestris DEGs and associated stats.
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Set variables based on argument.
  
  if (x == "brain") {
    controlColumns <- c(1, 24:31, 40)
    removalColumns <- c(1, 32:40)
  } else if (x == "fatbody") {
    controlColumns <- c(1, 24:31, 40)
    removalColumns <- c(1, 32:40)
  } else if (x == "ovaries") {
    controlColumns <- c(1, 26:33, 42)
    removalColumns <- c(1, 34:42)
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  # Load overall gene list based on argument.
  
  overallGeneList <- 
    read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Split overallGeneList into individual lists.
  
  # First, split by treatment, rename columns and remove NAs in expression_with_age.
  
  columnNames <- c("symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", 
                   "pvalue", "padj", "DEG", "expression_with_age", "name")
  
  controlResults <- overallGeneList[, controlColumns]
  
  colnames(controlResults) <- columnNames
  
  controlResults <- controlResults[!(is.na(controlResults[, "expression_with_age"])), ]
  
  removalResults <- overallGeneList[, removalColumns]
  
  colnames(removalResults) <- columnNames
  
  removalResults <- removalResults[!(is.na(removalResults[, "expression_with_age"])), ]
  
  # Split each treatment into up- and down-regulated DEGs.
  # Up-regulated = more expressed in older queens.
  # Down-regulated = more expressed in younger queens.
  
  controlUpDEGs <- 
    controlResults[controlResults[, "expression_with_age"] == "up-regulated", ]
  
  controlDownDEGs <- 
    controlResults[controlResults[, "expression_with_age"] == "down-regulated", ]
  
  removalUpDEGs <- 
    removalResults[removalResults[, "expression_with_age"] == "up-regulated", ]
  
  removalDownDEGs <- 
    removalResults[removalResults[, "expression_with_age"] == "down-regulated", ]
  
  # Remove DEG and expression_with_age columns.
  
  controlUpDEGs <- controlUpDEGs[, c(1:7, 10)]
  
  controlDownDEGs <- controlDownDEGs[, c(1:7, 10)]
  
  removalUpDEGs <- removalUpDEGs[, c(1:7, 10)]
  
  removalDownDEGs <- removalDownDEGs[, c(1:7, 10)]
  
  # Write results to .csv files.
  
  write.csv(controlUpDEGs, paste0(x, "_results_C_treatment_up_regulated_LFC0.csv"),
            row.names = FALSE)
  
  write.csv(controlDownDEGs, paste0(x, "_results_C_treatment_down_regulated_LFC0.csv"),
            row.names = FALSE)
  
  write.csv(removalUpDEGs, paste0(x, "_results_R_treatment_up_regulated_LFC0.csv"),
            row.names = FALSE)
  
  write.csv(removalDownDEGs, paste0(x, "_results_R_treatment_down_regulated_LFC0.csv"),
            row.names = FALSE)
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsBrain <- BuildDDSForDEAnalysis("brain")

ddsFatBody <- BuildDDSForDEAnalysis("fatbody")

ddsOvaries <- BuildDDSForDEAnalysis("ovaries")

# Make and Save Gene Lists ----

# Brain.

dir.create("../02_outputs/02_brain/31_DESeq2_DEG_lists")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../02_outputs/02_brain/31_DESeq2_DEG_lists")

# Make gene lists.

MakeGeneList("brain")

MakeDEGLists("brain")

# Fat body.

# Create directory.

dir.create("../../03_fatbody/31_DESeq2_DEG_lists")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists")

# Make gene lists.

MakeGeneList("fatbody")

MakeDEGLists("fatbody")

# Ovaries.

# Create directory.

dir.create("../../01_ovaries/31_DESeq2_DEG_lists")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists")

# Make gene lists.

MakeGeneList("ovaries")

MakeDEGLists("ovaries")
