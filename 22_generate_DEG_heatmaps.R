#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce a heatmap of differentially expressed genes for each
# tissue.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file, Kallisto abundances, virus_samples.csv,
# genesymbols2genenames.txt file.

# Outputs:
# .svg figure of a heatmap for each tissue.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(pheatmap)  # pheatmap().
library(ggplot2)  # ggsave().
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

DetermineTop50DEGs <- function (x) {
  # Generate data.frame of the 50 most highly differentially expressed genes (DEGs) 
  # (or all DEGs if fewer than 50 in total) for a given tissue.
  # 
  # Args:
  #   x: string denoting which tissue to generate a list of DEGs for ("brain",
  #      "fatbody" or "ovaries").
  #
  # Returns:
  #   A data.frame of the most highly differentially expressed genes.
  
  # Set variables based on argument.
  
  if (x == "brain") {
    dds <- ddsBrain
    reorderColumns <- c(12:22, 1:11)
    sampleColumns <- 8:29
  } else if (x == "fatbody") {
    dds <- ddsFatBody
    reorderColumns <- c(11:22, 1:10)
    sampleColumns <- 8:29
  } else if (x == "ovaries") {
    dds <- ddsOvaries
    reorderColumns <- c(13:24, 1:12)
    sampleColumns <- 8:31
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  # Generate data.frame of log transformed counts for each gene.
  
  # rlog transform the counts.
  
  ddsTransform <- rlog(dds)
  
  # Extract transformed counts into a data.frame.
  
  rlogNormCounts <- as.data.frame(assay(ddsTransform))
  
  # Reorder columns in rlogNormCounts.
  
  rlogNormCounts <-rlogNormCounts[, reorderColumns]
  
  # Add a column with the gene symbol to the data.frame.
  
  rlogNormCounts$symbol <- row.names(rlogNormCounts)
  
  # Extract the results of the DEG analysis for each gene, for both the control
  # and removal treatments.
  
  controlResults <- as.data.frame(results(dds, contrast = c("condition", "C_TP2G", "C_TP1G"),
                            alpha = 0.05, lfcThreshold = 0)) # Stats for control
  
  removalResults <- as.data.frame(results(dds, contrast = c("condition", "R_TP2G", "R_TP1G"),
                                          alpha = 0.05, lfcThreshold = 0)) # Stats for removal
  
  # Remove NAs from padj column (represent genes excluded from analysis as all 
  # counts were 0 or it contained an extreme count outlier).
  
  controlResults <- controlResults[!is.na(controlResults$padj), ]
  
  removalResults <- removalResults[!is.na(removalResults$padj), ]
  
  # Add gene symbol as a column.
  
  controlResults$symbol <- row.names(controlResults)
  
  removalResults$symbol <- row.names(removalResults)
  
  # Index only the significant results (padj < 0.05).
  
  sigControlResults <- controlResults[controlResults[, "padj"] < 0.05 , ]
  
  sigRemovalResults <- removalResults[removalResults[, "padj"] < 0.05 , ]
  
  # Combine Control and Removal results.
  
  allSigResults <- rbind(sigControlResults, sigRemovalResults)
  
  # Make all log2FoldChange values absolute, so that they can be ordered by
  # change irrespective of direction.
  
  allSigResults$log2FoldChange <- abs(allSigResults$log2FoldChange)
  
  # Sort data.frame so that rows are sorted by logFoldChange.
  
  sortedSigResults <- 
    allSigResults[order(-(allSigResults$log2FoldChange)), ]
  
  # Duplicated genes removed for data.frame, leaving the highest logFoldchange.
  
  nonDuplicateSigResults <- sortedSigResults[!duplicated(sortedSigResults$symbol), ]
  
  # Select top 50 DEGs by highest logFoldChange (or all DEGs if fewer).
  
  if (length(nonDuplicateSigResults$log2FoldChange) < 50) {
    numberOfGenes <- length(nonDuplicateSigResults$log2FoldChange)
  } else {
    numberOfGenes <- 50
  }
  
  # Index top DEGs.
  
  topDEGs <- nonDuplicateSigResults[1:numberOfGenes, ]
  
  # Combine topDEGs with the matrix of transformed counts
  
  plottingData <- merge(topDEGs, rlogNormCounts, all.x = TRUE)
  
  # Remove columns no longer needed.
  
  plottingData <- plottingData[, c(1, sampleColumns)]
  
  # Add gene names to symbols.
  
  plottingData <- merge(plottingData, symbolToName, all.x = TRUE)
  
  # Return plottingData.
  
  return(plottingData)
  
}

PlotHeatmap <- function (x) {
  # Takes a data.frame of expression data and generates a heatmap.
  #
  # Args:
  #   x: string denoting which tissue the data are for ("brain", "fatbody" 
  #      or "ovaries").
  #
  # Returns:
  #   A heatmap showing expression of the genes in the data.frame.
  
  # Set variables based on argument.
  
  if (x == "brain") {
    plottingData <- brainData
    y <- 2:23
    columnGaps <- c(6, 11, 17)
    virusNumbers <- c(13:23, 1:10, 12)
    timePointNumbers <- rep(c(rep("TP1G", times = 6), rep("TP2G", times = 5)), times = 2)
    treatmentNumbers <- c(rep("Removal queens", times = 11), rep("Control queens", times = 11))
  } else if (x == "fatbody") {
    plottingData <- fatBodyData
    y <- 2:23
    columnGaps <- c(6, 12, 18)
    virusNumbers <- c(13:24, 1:9, 11)
    timePointNumbers <- c(rep("TP1G", times = 6), rep("TP2G", times = 6), 
                          rep("TP1G", times = 6), rep("TP2G", times = 4))
    treatmentNumbers <- c(rep("Removal queens", times = 12), rep("Control queens", times = 10))
  } else if (x == "ovaries") {
    plottingData <- ovariesData
    y <- 2:25
    columnGaps <- c(6, 12, 18)
    virusNumbers <- c(13:24, 1:12)
    timePointNumbers <- rep(c(rep("TP1G", times = 6), rep("TP2G", times = 6)), times = 2)
    treatmentNumbers <- c(rep("Removal queens", times = 12), rep("Control queens", times = 12))
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  # Transform data.frame into a matrix.
  
  plotMatrix <- as.matrix(plottingData[, y])
  
  # Add gene names to the matrix.
  
  row.names(plotMatrix) <- plottingData[, "name"]
  
  # Generate keys and annotations.
  
  colAnnotations <- 
    data.frame(Virus = virus_samples[virusNumbers, 2],
               "Time-point" = timePointNumbers,
               Treatment = treatmentNumbers,
               check.names = FALSE)
  
  row.names(colAnnotations) <- colnames(plotMatrix)
  
  annotationColours <- 
    list(Virus = c("Yes" = "white", "No" = "black"),
         "Time-point" = c("TP1G" = "red", "TP2G" = "blue"),
         Treatment = c("Removal queens" = "orange", "Control queens" = "gray68"))
  
  # Plot heatmap.
  
  heatmapPlot <- pheatmap(plotMatrix,
                          scale = "row",
                          cluster_rows = TRUE,
                          cluster_cols = FALSE,
                          gaps_col = columnGaps,
                          annotation_col = colAnnotations,
                          annotation_colors = annotationColours,
                          legend = TRUE,
                          cellwidth=15,
                          cellheight=10,
                          fontsize = 7)
  
  # Return heatmap.
  
  return(heatmapPlot)
  
}

AddGeneSymbol <- function(x, y) {
  # Adds the gene symbol to the gene name for genes that share similar names.
  #
  # Args:
  #   x: data.frame containing a list of DEGs.
  #   y: number indicating the row where the gene symbol is to be added.
  #
  # Returns:
  #   The data.frame with the gene symbol added at the appropriate place.
  
  # Paste the gene symbol at the end of the gene name.
  
  x[y, "name"] <- paste0(x[y, "name"], " ", x[y, "symbol"])
  
  # Return x.
  
  return(x)
  
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

# Brain.

ddsBrain <- BuildDDSForDEAnalysis("brain")

# Fat body.

ddsFatBody <- BuildDDSForDEAnalysis("fatbody")

# Ovaries.

ddsOvaries <- BuildDDSForDEAnalysis("ovaries")

# Make Heatmaps ----

# Data for heatmaps.

# Brain.

brainData <- DetermineTop50DEGs("brain")

# Fat body.

fatBodyData <- DetermineTop50DEGs("fatbody")

# Ovaries.

ovariesData <- DetermineTop50DEGs("ovaries")

# Format gene names.
# Print gene names to check for duplicate names.
# Brain.

brainData$name

# Add gene symbols to genes with very similar names.

brainData <- AddGeneSymbol(brainData, 7)

brainData <- AddGeneSymbol(brainData, 8)

brainData <- AddGeneSymbol(brainData, 17)

brainData <- AddGeneSymbol(brainData, 22)

brainData <- AddGeneSymbol(brainData, 23)

brainData <- AddGeneSymbol(brainData, 26)

brainData <- AddGeneSymbol(brainData, 28)

brainData <- AddGeneSymbol(brainData, 32)

brainData <- AddGeneSymbol(brainData, 37)

brainData <- AddGeneSymbol(brainData, 50)

# Fat body.
# Print gene names to check for duplicate names.

fatBodyData$name

# Add gene symbols to genes with very similar names.

fatBodyData <- AddGeneSymbol(fatBodyData, 4)

fatBodyData <- AddGeneSymbol(fatBodyData, 17)

fatBodyData <- AddGeneSymbol(fatBodyData, 38)

fatBodyData <- AddGeneSymbol(fatBodyData, 44)

fatBodyData[48, 24] <- paste0("uncharacterized protein ", fatBodyData[48, 1])

# Ovaries.
# Print gene names to check for duplicate names.

ovariesData$name

# Add gene symbols to genes with very similar names.

ovariesData <- AddGeneSymbol(ovariesData, 6)

ovariesData <- AddGeneSymbol(ovariesData, 7)

ovariesData <- AddGeneSymbol(ovariesData, 28)

ovariesData[35, 26] <- "uncharacterized protein LOC100651091"

# Plot heatmaps.

# Brain.

brainHeatmapPlot <- PlotHeatmap("brain")

# Fat body.

fatBodyHeatmapPlot <- PlotHeatmap("fatbody")

# Ovaries.

ovariesHeatmapPlot <- PlotHeatmap("ovaries")

# Save Heatmaps ----

# Brain.

# Create directory.

dir.create("../02_outputs/02_brain/32_DEG_top50_heatmaps")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../02_outputs/02_brain/32_DEG_top50_heatmaps")

# Save plot.

ggsave("00_NER0008751_obj1_fig_S11_brain_heatmap.svg", brainHeatmapPlot,
       width = 29, height = 21, units = "cm")

# Fat body.

# Create directory.

dir.create("../../03_fatbody/32_DEG_top50_heatmaps")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../03_fatbody/32_DEG_top50_heatmaps")

# Save plot.

ggsave("00_NER0008751_obj1_fig_S12_fat_body_heatmap.svg", fatBodyHeatmapPlot,
       width = 29, height = 21, units = "cm")

# Ovaries.

# Create directory.

dir.create("../../01_ovaries/32_DEG_top50_heatmaps")
# Will produce a warning if directory already exists.

# Change directory.

setwd("../../01_ovaries/32_DEG_top50_heatmaps")

# Save plot.

ggsave("00_NER0008751_obj1_fig_S13_ovaries_heatmap.svg", ovariesHeatmapPlot,
       width = 29, height = 21, units = "cm")
