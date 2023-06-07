#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Plot read count vs. HISAT2 alignment for each tissue.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file, Kallisto pseudoalignment abundances, 
# and HISAT2 alignment stats.

# Outputs:
# .svg files of plots showing read count vs. overall sample alignment.
# .csv files of normalised read counts from Kallisto.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ggplot2)  # ggscatter() and associated functions.
library(ggpubr)  # ggarrange().
library(ggrepel)  # geom_text_repel().
# DESeq2 and tximport loaded via the scripts in the "Custom functions" section.

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Transcript to Gene File ----

setwd("../00_data/10_transcripts2genes")

t2g <- read.table(file = "Bter_v1_transcripts2genes.txt",
                  header = FALSE,
                  col.names = c("TXNAME",
                                "GENEID"))

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_SumT2GBuildDDS.R")

# FUNCTION DEFINITIONS ----

plotReadCountsVsAlignment <- function (x, y = "no") {
  # Plots read counts for Holobiont sequences against percentage of reads that 
  # align to Bombus terrestris genome with HISAT2.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("brain", "fatbody" 
  #      or "ovaries").
  #   y: whether plots for all 10 sequences should be returned ("no" - default),
  #      or only plots for slow bee paralysis virus (SBPV) ("yes") with name 
  #      labels for the samples. 
  #      
  # Returns:
  #   A ggplot object containing scatter plots.
  #
  # NOTE: The working directory needs to be the 02_outputs folder when using 
  #       this script.
  
  # Check working directory is correct.
  
  if (!(grepl("02_outputs", getwd()))) {
    stop('Working directory is incorrect, please set to "02_outputs".')
  }
  
  # Set arguments.
  
  if (x == "brain") {
    dds <- ddsBrain
    numberTissue <- "02_brain"
    outputFile <- "00_NER0008751_obj1_sup_file_Sz_brain_holobee_kallisto.csv"
    stringStart <- 13
    stringStop <- 23
  } else if (x == "fatbody") {
    dds <- ddsFatBody
    numberTissue <- "03_fatbody"
    outputFile <- "00_NER0008751_obj1_sup_file_Saa_fatbody_holobee_kallisto.csv"
    stringStart <- 15
    stringStop <- 25
  } else if (x == "ovaries") {
    dds <- ddsOvaries
    numberTissue <- "01_ovaries"
    outputFile <- "00_NER0008751_obj1_sup_file_Sab_ovary_holobee_kallisto.csv"
    stringStart <- 13
    stringStop <- 23
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  # Extract read counts into a data.frame.
  
  readCountsDF <- as.data.frame(assay(dds))
  
  # Output read count data.frame.
  
  write.csv(readCountsDF, 
            paste0(numberTissue, 
                  "/21_kallisto_pseudoalignment_holobee_summaries/", 
                   outputFile))
  
  # Transpose read count data.frame (swap rows and columns).
  
  readCountsDF2 <- data.frame(t(readCountsDF))
  
  # Add sample name column.
  
  readCountsDF2$Sample <- colnames(readCountsDF)
  
  # Import HISAT2 alignment stats for the specified tissue.
  
  setwd(paste0(numberTissue, "/02_quality_control_reports/", 
               "10_NER0008751_obj1_", x, "_hisat2_multiqc_data"))
  
  alignmentStats <- read.table("multiqc_hisat2.txt",
                                header = TRUE)
  
  # Reset working directory.
  
  setwd("../../../")
  
  # Add sample name to alignmentStats.
  
  alignmentStats$Sample <- substr(alignmentStats$Sample, 
                                  start = stringStart, 
                                  stop = stringStop)
  
  # Remove every other row, to reduce the duplication.
  
  alignmentStats <- 
    alignmentStats[(seq(from = 1, to = length(alignmentStats$Sample), by = 2)) ,]  
  
  # Merge alignment and read counts.
  
  mergedDF <- merge(readCountsDF2, alignmentStats)
  
  # Generate scatter plots for alignment vs. read count.
  
  scatter1 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.171673384.gb.EU569326.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter2 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.238800115.gb.GQ162109.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter3 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.425869910.gb.KC177809.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter4 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.312299480.gb.HQ619890.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter5 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.307148859.gb.HM237361.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter6 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.297578408.gb.GU938761.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter7 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.296005646.ref.NC_014137.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter8 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.374110828.gb.JN969316.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter9 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.374110809.gb.JN969337.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                         vjust = 0.5))
  
  scatter10 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                         y = "gi.374110808.gb.JN969338.1.",
                         xlab = "Overall sample alignment rate", 
                         ylab = "Normalised read counts for holobiont") +
                         theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                                          vjust = 0.5))
  
  scatter11 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.297578408.gb.GU938761.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5)) +
                        geom_text_repel(aes(label = Sample), size = 2, max.overlaps = Inf) +
                        theme(legend.position="none")
  
  scatter12 <- ggscatter(mergedDF, x = "overall_alignment_rate", 
                        y = "gi.296005646.ref.NC_014137.1.",
                        xlab = "Overall sample alignment rate", 
                        ylab = "Normalised read counts for holobiont") +
                        theme(axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5)) +
                        geom_text_repel(aes(label = Sample), size = 2, max.overlaps = Inf) +
                        theme(legend.position="none")
  
  # Combine plots into a single plot.
  
  if (y == "no") {
    # Combine all of the plots.
    
    scatterPlotsToReturn <- ggarrange(scatter1, scatter2, scatter3, scatter4, scatter5, 
                                      scatter6, scatter7, scatter8, scatter9, scatter10,
                                      labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"),
                                      ncol = 3, nrow = 4) 
  } else if (y == "yes") {
    # Only return the virus plots.
    
    scatterPlotsToReturn <- ggarrange(scatter11, scatter12, 
                                      ncol = 2, nrow = 1)
  }
  
  # Return plots.
  
  return(scatterPlotsToReturn)
  
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsBrain <- SumT2GBuildDDS("brain", "holobee")

ddsFatBody <- SumT2GBuildDDS("fatbody", "holobee")

ddsOvaries <- SumT2GBuildDDS("ovaries", "holobee")

# Generate Plots ----

# Set working directory.

setwd("../02_outputs/")

# All plots.

brainPlots <- plotReadCountsVsAlignment ("brain")

fatBodyPlots <- plotReadCountsVsAlignment ("fatbody")

ovariesPlots <- plotReadCountsVsAlignment ("ovaries")

# Just virus plots for each tissue.

brainVirusPlots <- plotReadCountsVsAlignment ("brain", "yes")

fatBodyVirusPlots <- plotReadCountsVsAlignment ("fatbody", "yes")

ovariesVirusPlots <- plotReadCountsVsAlignment ("ovaries", "yes")

allVirusPlots <- ggarrange(brainVirusPlots,
                            fatBodyVirusPlots,
                            ovariesVirusPlots,
                            ncol = 1, nrow = 3)

# Save Plots ----

# Set working directory.

setwd("00_all_tissues/")

ggsave("10_NER0008751_obj1_fig_S16_fat_body_holobee_kallisto.svg", 
       fatBodyPlots,
       height = 29, width = 21, units = "cm")

ggsave("11_NER0008751_obj1_fig_S17_virus_all_tissues.svg", 
       allVirusPlots,
       height = 29, width = 21, units = "cm")              
