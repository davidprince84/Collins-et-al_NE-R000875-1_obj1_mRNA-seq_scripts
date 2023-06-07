#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Combine gene ontology (GO) enrichment results into a single table.
#-------------------------------------------------------------------------------
# Inputs:
# Results of the GO enrichment results from script #23.

# Outputs:
# A .txt file combining the results for enriched Biological Processes GO terms.
# NOTE: .txt file used rather than .csv file to prevent Excel from automatically
# formatting GeneRatio and BgRatio columns, which can result in dates being 
# inserted rather than fractions.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# No packages loaded.

# LOAD DATA ----

# Data loaded via function.

# FUNCTION DEFINITIONS ----

LoadGOResults <- function (x, y, z) {
  # Load the results of a specified GO enrichment analysis and formats them so
  # that they can be combined into a single table.
  #
  # Args:
  #   x: string denoting which tissue to load enriched GO term list from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting which treatment the B. terrestris enriched GO term list
  #      is from ("C" or "R").
  #   z: string denoting whether the enriched GO term list to load is for up-regulated 
  #      genes (genes more highly expressed in time point 2 compared to time 
  #      point one) or down-regulated genes (genes more highly expressed in time
  #      point 1 compared to time point 2) ("up" or "down").
  #
  # Returns:
  #   A data.frame containing the formatted results.
  #
  # Note: The working directory needs to be set to 
  # 40_GO_enrichment_analysis for the relevant tissue.
  
  # Load data based on arguments.
  
  goResults <- 
    read.table(paste(x, y, z, 
                     "regulated_DEGs_enriched_GO_results.txt", sep = "_"),
               header = TRUE)
  
  # Set column names.
  
  columnNames <- c("GO term ID", "GO term description", "Gene ratio",
                   "Background gene ratio", "p value", "Adjusted p value",
                   "q value", "Flybase gene ID", "Count", "Redundant")
  
  columnOrder <- c("Tissue", "Treatment", "Expression with age", "Ontology",
                   columnNames)
  
  # Check that enriched GO terms were found and act accordingly.
  
  if (goResults[1,1] == "No enriched GO terms found.") {
    
    formattedGOResults <- data.frame(Tissue = x,
                                     Treatment = y,
                                     "Expression with age" = paste0(z, "-regulated"),
                                     Ontology = "Biological Process",
                                     "GO term ID" = "NA",
                                     "GO term description" = "No enriched GO terms found",
                                     "Gene ratio" = "NA",
                                     "Background gene ratio" = "NA",
                                     "p value" = "NA",
                                     "Adjusted p value" = "NA",
                                     "q value" = "NA",
                                     "Flybase gene ID" = "NA",
                                     "Redundant" = "NA",
                                     check.names = FALSE)
  } else {
    
    # Rename columns.
    
    colnames(goResults) <- columnNames
    
    # Add additional columns to the data.
    
    goResults$Tissue <- rep(x, times = length(goResults$`GO term ID`))
    
    goResults$Treatment <- rep(y, times = length(goResults$`GO term ID`))
    
    goResults$"Expression with age" <- rep(paste0(z, "-regulated"), times = length(goResults$`GO term ID`))
    
    goResults$Ontology <- rep("Biological Process", times = length(goResults$`GO term ID`))
    
    # Reorder columns.
    
    formattedGOResults <- goResults[, columnOrder[c(1:12, 14)]]  # Count column removed.
    
  }
  
  # Return formattedGOResults.
  
  return(formattedGOResults)
  
}

# EXECUTED STATEMENTS ----

# Load and Format Brain GO Results ----

# Set working directory.

setwd("../02_outputs/02_brain/40_GO_enrichment_analysis/")

# Load and format results.

# C results.

brainCUpResults <- LoadGOResults("brain", "C", "up")

brainCDownResults <- LoadGOResults("brain", "C", "down")

# R results.

brainRUpResults <- LoadGOResults("brain", "R", "up")

brainRDownResults <- LoadGOResults("brain", "R", "down")

# Load and Format Fat Body GO Results ----

# Set working directory.

setwd("../../03_fatbody/40_GO_enrichment_analysis/")

# Load and format results.

# C results.

fatBodyCUpResults <- LoadGOResults("fatbody", "C", "up")

fatBodyCDownResults <- LoadGOResults("fatbody", "C", "down")

# R results.

fatBodyRUpResults <- LoadGOResults("fatbody", "R", "up")

fatBodyRDownResults <- LoadGOResults("fatbody", "R", "down")

# Load and Format Ovaries GO Results ----

# Set working directory.

setwd("../../01_ovaries/40_GO_enrichment_analysis/")

# Load and format results.

# C results.

ovariesCUpResults <- LoadGOResults("ovaries", "C", "up")

ovariesCDownResults <- LoadGOResults("ovaries", "C", "down")

# R results.

ovariesRUpResults <- LoadGOResults("ovaries", "R", "up")

ovariesRDownResults <- LoadGOResults("ovaries", "R", "down")

# Combine All Results into a Single Data.Frame ----

combinedResults <- rbind(brainRUpResults, brainRDownResults,
                         brainCUpResults, brainCDownResults,
                         fatBodyRUpResults, fatBodyRDownResults,
                         fatBodyCUpResults, fatBodyCDownResults, 
                         ovariesRUpResults, ovariesRDownResults,
                         ovariesCUpResults, ovariesCDownResults)

# Save the Combined Results ----

# Set working directory.

setwd("../../00_all_tissues/")

# Save as a .txt file.

write.table(combinedResults, 
            "21_NER0008751_obj1_table_S9_GO_enrichment_results.txt", 
            sep = "\t", row.names = FALSE)
