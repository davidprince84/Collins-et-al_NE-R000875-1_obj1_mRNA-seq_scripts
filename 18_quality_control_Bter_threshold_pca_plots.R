#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Produce principal component analysis and normalisation boxplots from 
# Kallisto data with technical replicates collapsed, read number threshold 
# added, and presence of virus-aligning reads included or excluded in the model.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file, Kallisto abundances, virus_samples.csv

# Outputs:
# .svg figures of the plots.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(reshape2)  # melt()
library(ggpubr)  # ggarrange(), annotate_figure()
# ggrepel, ggplot2, DESeq2 and tximport loaded via the scripts in the 
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

# Virus Samples ----

setwd("../11_virus_samples")

virus_samples <- read.csv(file = "virus_samples.csv",
                          col.names = c("sample", "virus"))

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_BuildDDSForDEAnalysis.R")

source("f_CustomPCAPlot.R")

# FUNCTION DEFINITIONS ----

GenerateExploratoryPlots <- function(x, y = "no") {
  # Generates a figure showing the normalisation of the gene expression data by 
  # DESeq2 as a boxplot and the clustering of this data as a principal component
  # analysis (PCA) plot.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("brain", "fatbody" 
  #      or "ovaries").
  #   y: string denoting whether samples containing many virus-aligning
  #      reads should be excluded from the analysis ("yes" and "no" (default)).
  #
  # Returns:
  #   A ggplot annotated figure with 2 subplots arranged by the ggpubr package. 
  
  # Set variables based on argument.
  
  if (x == "brain") {
    ddsName <- "ddsBrain"
    pcaName <- "pcaBrain"
    tissue <- "Brain"
  } else if (x == "fatbody") {
    ddsName <- "ddsFatBody"
    pcaName <- "pcaFatBody"
    tissue <- "Fat body"
  } else if (x == "ovaries") {
    ddsName <- "ddsOvaries"
    pcaName <- "pcaOvaries"
    tissue <- "Ovaries"
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  if (y == "no") {
    samplesName <- "AllSamples"
  } else if (y == "yes") {
    samplesName <- "NoVirus"
  } else {
    stop('Argument y must be "yes" or "no" (default).')
  }
  
  dds <- eval(parse(text = paste0(ddsName, samplesName)))
  
  pca <- eval(parse(text = paste0(pcaName, samplesName)))
  
  # Transform the data to stabilize the variance across the mean.
  # Select  rlog transformation chosen if there are relatively few 
  # samples (n < 30), or vst transformation if there are more than 30
  # samples.
  
  if (length(dds$sample) < 30) {
    ddsTransform <- rlog(dds)
    yLabelContents <- "rlog transformed counts"
  } else if (length(dds$sample) >= 30) {
    ddsTransform <- vst(dds)
    yLabelContents <- "vst transformed counts"
  }
  
  # Box plot of normalized, transformed counts.
  
  # Format data for plot.
  
  boxPlotData <- melt(assay(ddsTransform))
  colnames(boxPlotData) <- c("gene", "sample", "transformed_counts")
  
  # Calculate median value of counts.
  
  samplesMedian <- median(boxPlotData$transformed_counts)
  
  # Generate box plot.
  
  boxPlot <- ggplot(boxPlotData, aes(x = sample, y = transformed_counts)) + 
    geom_boxplot() +
    geom_hline(yintercept = samplesMedian) +  # Add median line
    ylab(yLabelContents)        
  
  # Format box plot.
  
  formattedBoxPlot <- boxPlot + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, 
                                     vjust = 0.5))
  
  # Group plots into a single plot and annotate with tissue.
  
  figure <- ggarrange(formattedBoxPlot, pca, 
                      labels = c("a", "b"),
                      ncol = 2, nrow = 1)
  
  figure <- annotate_figure(figure,
                            top = tissue)
  
  # Return figure.
  
  return(figure)
  
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue With All Samples ----

ddsBrainAllSamples <- BuildDDSForDEAnalysis("brain")

ddsFatBodyAllSamples <- BuildDDSForDEAnalysis("fatbody")

ddsOvariesAllSamples <- BuildDDSForDEAnalysis("ovaries")

# Build dds Objects for Each Tissue With "With Virus" Samples Excluded ----

ddsBrainNoVirus <- BuildDDSForDEAnalysis("brain", y = "no_virus")

ddsFatBodyNoVirus <- BuildDDSForDEAnalysis("fatbody", y = "no_virus")

ddsOvariesNoVirus <- BuildDDSForDEAnalysis("ovaries", y = "no_virus")

# PCA Calls ----

pcaBrainAllSamples <- CustomPCAPlot("brain")

pcaFatBodyAllSamples <- CustomPCAPlot("fatbody")

pcaOvariesAllSamples <- CustomPCAPlot("ovaries")

pcaBrainNoVirus <- CustomPCAPlot("brain", y = "yes")

pcaFatBodyNoVirus <- CustomPCAPlot("fatbody", y = "yes")

pcaOvariesNoVirus <- CustomPCAPlot("ovaries", y = "yes")

# Make Figures With Normalisation Boxplots ----

brainPlotsAllSamples <- GenerateExploratoryPlots("brain")

fatBodyPlotsAllSamples <- GenerateExploratoryPlots("fatbody")

ovariesPlotsAllSamples <- GenerateExploratoryPlots("ovaries")

brainPlotsNoVirus <- GenerateExploratoryPlots("brain", "yes")

fatBodyPlotsNoVirus <- GenerateExploratoryPlots("fatbody", "yes")

ovariesPlotsNoVirus <- GenerateExploratoryPlots("ovaries", "yes")

# Save Plots ----

# Brain.

setwd("../02_outputs/02_brain/02_quality_control_reports")

ggsave("21_NER0008751_obj1_fig_S20_brain_exploratory_plots_all_samples.svg",
       brainPlotsAllSamples, width = 29.7, height = 21, units = "cm")

ggsave("22_NER0008751_obj1_brain_exploratory_plots_no_virus.svg",
       brainPlotsNoVirus, width = 29.7, height = 21, units = "cm")

# Fat body.

setwd("../../03_fatbody/02_quality_control_reports")

ggsave("21_NER0008751_obj1_fig_S21_fat_body_exploratory_plots_all_samples.svg",
       fatBodyPlotsAllSamples, width = 29.7, height = 21, units = "cm")

ggsave("22_NER0008751_obj1_fat_body_exploratory_plots_no_virus.svg",
       fatBodyPlotsNoVirus, width = 29.7, height = 21, units = "cm")

# Ovaries.

setwd("../../01_ovaries/02_quality_control_reports")

ggsave("21_NER0008751_obj1_fig_S22_ovaries_exploratory_plots_all_samples.svg",
       ovariesPlotsAllSamples, width = 29.7, height = 21, units = "cm")

ggsave("22_NER0008751_obj1_fig_ovaries_exploratory_plots_no_virus.svg",
       ovariesPlotsNoVirus, width = 29.7, height = 21, units = "cm")
