#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 17.02.2021
# File last updated: 03.01.2023
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq and comparative analysis.
# Tasks: Identify single-copy orthologues from OrthoFinder results.
#-------------------------------------------------------------------------------
# Inputs:
# OrthoFinder results.

# Outputs:
# Data.frame of single-copy orthologues between the designated species.
#-------------------------------------------------------------------------------

# Note: This function requires the working directory to be set to the Orthologues
# subdirectory of the OrthoFinder results.

# FUNCTION DEFINITION ----

IdentifySingleCopyOrthos <- function (x, y) {
  # Identifies single copy orthologues from the OrthoFinder results for 
  # the two specified species in the arguments.
  #
  # Args:
  #   x: string denoting the first species name ("Apis_mellifera", 
  #      Bombus_terrestris" or "Drosophila_melanogaster").
  #   y: string denoting the second species name ("Apis_mellifera", 
  #      Bombus_terrestris" or "Drosophila_melanogaster").
  #
  # Returns:
  #   A data.frame of single-copy orthologues between the designated species.
  
  # Check that the same species name hasn't been given twice.
  
  if (x == y) {
    stop('Arguments x and y are the same. They must be different.')
  }
  
  # Check that the working directory in the OrthoFinder results.
  
  if (!(grepl("01_results/Orthologues", getwd()))) {
    stop('Error - please change working directory to "02_outputs/10_orthofinder/01_results/Orthologues".')
  }
  
  # Set working directory.
  
  setwd(paste0("Orthologues_", x))
  
  # Load the results specified by the arguments.
  
  orthoList <- 
    read.table(paste0(x, "__v__", y, ".tsv"), 
               sep = '\t', header = TRUE)
  
  # Identify single copy orthologuess between species x and y.
  
  # Remove rows where the species x column has multiple entries
  # (by searching for the comma that separates the entries).
  
  orthoListMinusXMulti <- 
    orthoList[!(grepl(", ", orthoList[, x], fixed = TRUE)), ]
  
  # Remove rows where the species y column has multiple entries
  # (by searching for the comma that separates the entries).
  
  singleCopyOrthoList <- 
    orthoListMinusXMulti[!(grepl(", ", orthoListMinusXMulti[, y], fixed = TRUE)), ]
  
  # Reset working directory.
  
  setwd("../")
  
  # Return singleCopyOrthoList.
  
  return(singleCopyOrthoList)
  
}
