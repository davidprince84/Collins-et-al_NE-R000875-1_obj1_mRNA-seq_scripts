#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Produce principal component analysis of technical replicates from 
# Kallisto data.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file and Kallisto abundances.

# Outputs:
# .svg figures of the plots.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

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

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_SumT2GBuildDDS.R")

source("f_CustomPCAPlot.R")

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsBrainAllSamples <- SumT2GBuildDDS(x = "brain", y = "bter", z = "no")

ddsFatBodyAllSamples <- SumT2GBuildDDS(x = "fatbody", y = "bter", z = "no")

ddsOvariesAllSamples <- SumT2GBuildDDS(x = "ovaries", y = "bter", z = "no")

# PCA Calls ----

pcaBrain <- CustomPCAPlot("brain", techreps = "uncollapsed")

pcaFatBody <- CustomPCAPlot("fatbody", techreps = "uncollapsed")

pcaOvaries <- CustomPCAPlot("ovaries", techreps = "uncollapsed")

# Save Graphs ----
# NOTE: The graphs are saved as separate PDFs, as combining the three different
# graphs onto one page makes it difficult to see which labels are
# attached to which points.

# Ovaries.

setwd("../02_outputs/01_ovaries/02_quality_control_reports")

ggsave("20_NER0008751_obj1_ovaries_deseq2_tech_rep_pca.pdf", 
       pcaOvaries, width = 29.7, height = 21, units = "cm")

# Brain.

setwd("../../02_brain/02_quality_control_reports")

ggsave("20_NER0008751_obj1_brain_deseq2_tech_rep_pca.pdf", 
       pcaBrain, width = 29.7, height = 21, units = "cm")

# Fat body.

setwd("../../03_fatbody/02_quality_control_reports")

ggsave("20_NER0008751_obj1_fatbody_deseq2_tech_rep_pca.pdf", 
       pcaFatBody, width = 29.7, height = 21, units = "cm")

# Remove dds objects to avoid confusion with other scripts.

rm(ddsBrainAllSamples, ddsFatBodyAllSamples, ddsOvariesAllSamples)
