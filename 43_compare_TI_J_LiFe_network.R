#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Compare Bombus terrestris (Bter) differentially 
# expressed genes (DEGs) with Drosophila melanogaster (Dmel) TI-J-LiFe network 
# genes, as defined by Korb et al. (2021).
#-------------------------------------------------------------------------------
# Inputs:
# Bter DEGs from the current study. OrthoFinder results. List of Dmel
# TI-J-LiFe genes from Korb et al. (2021).

# Outputs:
# .csv file recording all the results of the statistical tests.
# .csv file recording all the genes that overlap between the TI-J-LiFe genes and 
# the DEGs from the current study.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD CUSTOM FUNCTIONS ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

source("f_IdentifySingleCopyOrthos.R")

source("f_CompareBterDEGsAndDmelGenes.R")

# LOAD DATA ----

# Load Korb et al. (2021) Data ----

# Set working directory.

setwd("../00_data/06_gene_list_csv")

# Load data.

TIJLiFeGenes <- read.csv("korb_2021_table_s1.csv",
                         stringsAsFactors = FALSE)

# Set working directory to OrthoFinder results.

setwd("../../02_outputs/10_orthofinder/01_results/Orthologues/")

# Identify single-copy orthologues between Dmel and Bter.

dmelBterOrthologues <- IdentifySingleCopyOrthos("Drosophila_melanogaster",
                                                "Bombus_terrestris")

# Rename columns of OrthoFinder results to allow merge.

colnames(dmelBterOrthologues) <- 
  c("Orthogroup", "Flybase_gene_ID", "Bter_gene_ID")

# FUNCTION DEFINITIONS ----

# No functions defined.

# EXECUTED STATEMENTS ----

# Brain Comparison ----

# Set working directory.

setwd("../../../02_brain/31_DESeq2_DEG_lists/")

# Compare gene lists.

brainControlAllResults <- CompareBterDEGsAndDmelGenes("brain", "C", a = "life")

brainControlTop50Results <- CompareBterDEGsAndDmelGenes("brain", "C", 50, a = "life")

brainControlTop100Results <- CompareBterDEGsAndDmelGenes("brain", "C", 100, a = "life")

brainControlTop200Results <- CompareBterDEGsAndDmelGenes("brain", "C", 200, a = "life")

brainControlTop300Results <- CompareBterDEGsAndDmelGenes("brain", "C", 300, a = "life")
# Not enough DEGs to calculate "top500".

brainRemovalAllResults <- CompareBterDEGsAndDmelGenes("brain", "R", a = "life")
# Not enough DEGs to calculate "top50" or higher.

# Fat Body Comparison ----

# Set working directory.

setwd("../../03_fatbody/31_DESeq2_DEG_lists/")

# Compare gene lists.

fatBodyControlAllResults <- CompareBterDEGsAndDmelGenes("fatbody", "C", a = "life")

fatBodyControlTop50Results <- CompareBterDEGsAndDmelGenes("fatbody", "C", 50, a = "life")

fatBodyControlTop100Results <- CompareBterDEGsAndDmelGenes("fatbody", "C", 100, a = "life")

fatBodyControlTop200Results <- CompareBterDEGsAndDmelGenes("fatbody", "C", 200, a = "life")

fatBodyControlTop300Results <- CompareBterDEGsAndDmelGenes("fatbody", "C", 300, a = "life")

fatBodyControlTop500Results <- CompareBterDEGsAndDmelGenes("fatbody", "C", 500, a = "life")

fatBodyRemovalAllResults <- CompareBterDEGsAndDmelGenes("fatbody", "R", a = "life")

fatBodyRemovalTop50Results <- CompareBterDEGsAndDmelGenes("fatbody", "R", 50, a = "life")

fatBodyRemovalTop100Results <- CompareBterDEGsAndDmelGenes("fatbody", "R", 100, a = "life")

fatBodyRemovalTop200Results <- CompareBterDEGsAndDmelGenes("fatbody", "R", 200, a = "life")

fatBodyRemovalTop300Results <- CompareBterDEGsAndDmelGenes("fatbody", "R", 300, a = "life")
# Not enough DEGs to calculate "top500".

# Ovaries Comparison ----

# Set working directory.

setwd("../../01_ovaries/31_DESeq2_DEG_lists/")

# Compare gene lists.

ovariesControlAllResults <- CompareBterDEGsAndDmelGenes("ovaries", "C", a = "life")

ovariesControlTop50Results <- CompareBterDEGsAndDmelGenes("ovaries", "C", 50, a = "life")

ovariesControlTop100Results <- CompareBterDEGsAndDmelGenes("ovaries", "C", 100, a = "life")

ovariesControlTop200Results <- CompareBterDEGsAndDmelGenes("ovaries", "C", 200, a = "life")

ovariesControlTop300Results <- CompareBterDEGsAndDmelGenes("ovaries", "C", 300, a = "life")

ovariesControlTop500Results <- CompareBterDEGsAndDmelGenes("ovaries", "C", 500, a = "life")

ovariesRemovalAllResults <- CompareBterDEGsAndDmelGenes("ovaries", "R", a = "life")
# Not enough DEGs to calculate "top50" or higher.

# Combine Results and Save ----

# Combine results for stats.

combinedStats <- rbind(brainControlTop50Results$stats, brainControlTop100Results$stats,
                       brainControlTop200Results$stats, brainControlTop300Results$stats,
                       brainControlAllResults$stats, brainRemovalAllResults$stats,
                       fatBodyControlTop50Results$stats, fatBodyControlTop100Results$stats,
                       fatBodyControlTop200Results$stats, fatBodyControlTop300Results$stats,
                       fatBodyControlTop500Results$stats, fatBodyControlAllResults$stats,
                       fatBodyRemovalTop50Results$stats, fatBodyRemovalTop100Results$stats, 
                       fatBodyRemovalTop200Results$stats, fatBodyRemovalTop300Results$stats,
                       fatBodyRemovalAllResults$stats,
                       ovariesControlTop50Results$stats, ovariesControlTop100Results$stats,
                       ovariesControlTop200Results$stats, ovariesControlTop300Results$stats,
                       ovariesControlTop500Results$stats, ovariesControlAllResults$stats,
                       ovariesRemovalAllResults$stats)

# Combine results for overlapping genes.

combinedGenes <- cbind(brainControlTop50Results$genes, brainControlTop100Results$genes,
                       brainControlTop200Results$genes, brainControlTop300Results$genes,
                       brainControlAllResults$genes, brainRemovalAllResults$genes,
                       fatBodyControlTop50Results$genes, fatBodyControlTop100Results$genes,
                       fatBodyControlTop200Results$genes, fatBodyControlTop300Results$genes,
                       fatBodyControlTop500Results$genes, fatBodyControlAllResults$genes,
                       fatBodyRemovalTop50Results$genes, fatBodyRemovalTop100Results$genes, 
                       fatBodyRemovalTop200Results$genes, fatBodyRemovalTop300Results$genes,
                       fatBodyRemovalAllResults$genes,
                       ovariesControlTop50Results$genes, ovariesControlTop100Results$genes,
                       ovariesControlTop200Results$genes, ovariesControlTop300Results$genes,
                       ovariesControlTop500Results$genes, ovariesControlAllResults$genes,
                       ovariesRemovalAllResults$genes)

# Save results.

setwd("../../00_all_tissues/")

write.csv(combinedStats, "33_NER0008751_obj1_table_S14_TI-J-LiFe_stats_results.csv",
          row.names = FALSE)

write.csv(combinedGenes, "34_NER0008751_obj1_table_S15_TI-J-LiFe_overlapping_genes.csv",
          row.names = FALSE)
