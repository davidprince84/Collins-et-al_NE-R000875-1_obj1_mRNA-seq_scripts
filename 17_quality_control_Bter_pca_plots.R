#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Produce principal component analysis from Kallisto data with 
# technical replicates collapsed and symbols denoting virus-containing 
# samples.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file, Kallisto abundances, virus_samples.csv

# Outputs:
# .svg figure of the plots.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

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

virus_samples <- read.csv(file = "virus_samples.csv")

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_SumT2GBuildDDS.R")

source("f_CustomPCAPlot.R")

# FUNCTION DEFINITIONS ----

# No function defined.

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsBrainAllSamples <- SumT2GBuildDDS("brain")

ddsFatBodyAllSamples <- SumT2GBuildDDS("fatbody")

ddsOvariesAllSamples <- SumT2GBuildDDS("ovaries")

# PCA Calls ----

pcaBrain <- CustomPCAPlot("brain", virussymbols = "yes")

pcaFatBody <- CustomPCAPlot("fatbody", virussymbols = "yes")

pcaOvaries <- CustomPCAPlot("ovaries", virussymbols = "yes")

# Combine Plots to Make Figure ----

pcaFigure <- ggarrange(pcaBrain, pcaFatBody, pcaOvaries,
                       labels = c("a        Brain",
                                  "b        Fat body",
                                  "c        Ovaries"),
                       ncol = 2, nrow = 2)

# Save Figure ----

setwd("../02_outputs/00_all_tissues")

ggsave("12_NER0008751_obj1_fig_S18_virus_PCA.svg", pcaFigure, width = 29.7, 
       height = 21, units = "cm")

# Remove dds objects to avoid confusion with other scripts.

rm(ddsBrainAllSamples, ddsFatBodyAllSamples, ddsOvariesAllSamples)
