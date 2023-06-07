#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Prepare data for comparison of Bombus terrestris (Bter) differentially 
# expressed genes (DEGs) with Drosophila melanogaster (Dmel) GenAge data set.
#-------------------------------------------------------------------------------
# Inputs:
# GenAge Dmel gene list (from script #07). OrthoFinder results. Lists of Dmel
# DEGs from Chen et al. (2014) and Pacifico et al. (2018).

# Outputs:
# Data.frame object of GenAge genes with Bter orthologues ready for comparative 
# analysis. Data.frames of DEGs from Chen et al. (2014) and Pacifico et al. 
# (2018) with FlyBase IDs added.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(org.Dm.eg.db)  # mapIds()

# LOAD CUSTOM FUNCTIONS ----

source("f_IdentifySingleCopyOrthos.R")

# LOAD DATA ----

# Load GenAge Gene List ----

# Set working directory.
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

setwd("../00_data/06_gene_list_csv/")

# Load the gene list.

genAgeGenes <- read.csv("Dmel_genage.csv")

# Load Gene Lists from Dmel mRNA-seq Studies ----

# Chen et al. (2014).

# Load up-regulated DEGs.

chenUpDEGs <- read.csv("chen_2014_table_S1.csv", stringsAsFactors = FALSE,
                        skip = 1, nrows = 1866)

# Load down-regulated DEGs.

chenDownDEGs <- read.csv("chen_2014_table_S1.csv", stringsAsFactors = FALSE,
                          skip = 1870)

# Pacifico et al. (2018).

# 30d vs 5d data are broadly equivalent to the current study.

pacificoData <- read.csv("pacifico_2018_table_S1_female.csv",
                         skip = 3,
                         stringsAsFactors = FALSE)

# FUNCTION DEFINITIONS ----

ReformatAndAddFlybaseIDs <- function (x) {
  # Reformats the data.frame to remove unneccesary columns, adds the matching 
  # Flybase gene ID, and removes genes with no FlyBase ID.
  #
  # Args:
  #   x: a data.frame containing a column named "Dmel_gene_symbol".
  #
  # Returns:
  #   A reformatted data.frame with FlyBase IDs added to the genes.
  
  # Set variable based on x argument.
  
  if (grepl("chen", deparse(substitute(x)))) {
    columnsToRetain <- c(5, 6, 8)
  } else if (grepl("pacifico", deparse(substitute(x)))) {
    columnsToRetain <- c(1, 7, 8)
  } else {
    stop('Argument x must contain "chen" or "pacifico.')
  }
  
  # Retain relevant columns.
  
  fewerColumnDF <- x[, columnsToRetain]
  
  # Rename columns.
  # LFC = log fold-change, FDR = false discovery rate.
  
  colnames(fewerColumnDF) <- c("Dmel_gene_symbol", "LFC", "FDR")
  
  # Map Flybase gene ID by gene symbols.
  
  fewerColumnDF$Dmel_gene_ID <- mapIds(org.Dm.eg.db,
                                      keys = fewerColumnDF$Dmel_gene_symbol,
                                      column = "FLYBASE",
                                      keytype = "ALIAS")
  
  # Remove genes that do not have a FlyBase ID.
  
  retainedGenes <- fewerColumnDF[!(is.na(fewerColumnDF$Dmel_gene_ID)), ]
  
  # Return retainedGenes.
  
  return(retainedGenes)
  
}

# EXECUTED STATEMENTS ----

# Reformat GenAge Data for Comparisons ----

# Remove genes for which "longevity.influence" is "Unclear" or "Unannotated",
# and which have no entrez.gene.id.

genAgeGenes <- genAgeGenes[!(genAgeGenes$longevity.influence == "Unclear"), ]  # 3 genes.

genAgeGenes <- genAgeGenes[!(genAgeGenes$longevity.influence == "Unannotated"), ]  # 1 gene.

genAgeGenes <- genAgeGenes[!(is.na(genAgeGenes$entrez.gene.id)), ]  # 1 gene.

# Add Flybase gene ids to entrez gene ids of GenAge genes in order to compare
# list with Orthofinder results.

# Need to convert entrez.gene.id into a character vector.

genAgeGenes$entrez.gene.id <- as.character(genAgeGenes$entrez.gene.id)

# Map Flybase gene ID by entrez gene IDs.
# Column becomes a list if as.character() not used.
# But don't use the symbol rather than entrez.gene.id, as some of the symbols 
# are replicated and give the wrong gene ID.

genAgeGenes$Flybase_gene_ID <- as.character(mapIds(org.Dm.eg.db,
                                            keys = genAgeGenes$entrez.gene.id,
                                            column = "FLYBASE",
                                            keytype = "ENTREZID"))

# Remove genes with no Flybase gene ID.

genAgeGenes <- genAgeGenes[!(is.na(genAgeGenes$Flybase_gene_ID)), ]  # 3 genes.

# Reformat OrthoFinder Results ----

# Set working directory to OrthoFinder results.

setwd("../../02_outputs/10_orthofinder/01_results/Orthologues/")

# Identify single-copy orthologues between Dmel and Bter.

dmelBterOrthologues <- IdentifySingleCopyOrthos("Drosophila_melanogaster",
                                                "Bombus_terrestris")

# Combine GenAge Genes with OrthoFinder Results ----

# Rename columns of OrthoFinder results to allow merge.

colnames(dmelBterOrthologues) <- 
  c("Orthogroup", "Flybase_gene_ID", "Bter_gene_ID")

# Merge the two data.frames.

genAgeOrtho <- merge(genAgeGenes, dmelBterOrthologues)
# Results in 124 Dmel GenAge genes with Bter single-copy orthologues.

# Reformat Chen et al. (2014) Data ----

# Up-regulated DEGs.
# Up-regulated = more expressed at the second time point compared to the first,
# hence up-regulated with age.

chenUpDEGsFB <- ReformatAndAddFlybaseIDs(chenUpDEGs)

# Down-regulated DEGs.
# Down-regulated = more expressed at the first time point compared to the second,
# hence down-regulated with age.

chenDownDEGsFB <- ReformatAndAddFlybaseIDs(chenDownDEGs)

# Reformat Pacifico et al. (2018) Data ----

# Reformat all DEGs.

pacificoDataFB <- ReformatAndAddFlybaseIDs(pacificoData)

# Remove rows where FDR = NA.
# These are reported by DESeq2 if no test was applied (all counts for gene were
# 0), or the gene was excluded from the analysis because it was an extreme 
# count outlier.

pacificoDataFB <- pacificoDataFB[!is.na(pacificoDataFB$FDR), ]

# Remove rows where FDR < 0.05 (i.e. non-significant).
# These rows are possible because three time step comparisons were 
# presented together in the spreadsheet, but only one comparison is 
# relevant for the current study.

pacificoDEGsFB <- pacificoDataFB[pacificoDataFB[, "FDR"] < 0.05, ]

# Index the data.frame to separate up-regulated and down-regulated DEGs.
# Up-regulated = more expressed at the second time point compared to the first,
# hence up-regulated with age.
# Down-regulated = more expressed at the first time point compared to the second,
# hence down-regulated with age.

pacificoUpDEGsFB <- pacificoDEGsFB[pacificoDEGsFB[, "LFC"] > 0, ]

pacificoDownDEGsFB <- pacificoDEGsFB[pacificoDEGsFB[, "LFC"] < 0, ]

# Remove Objects No Longer Required ----

rm(genAgeGenes, IdentifySingleCopyOrthos,chenUpDEGs, 
   chenDownDEGs, pacificoData, pacificoDataFB, pacificoDEGsFB,
   ReformatAndAddFlybaseIDs)

# Reset Working Directory ----

setwd("../../../../01_scripts/")
