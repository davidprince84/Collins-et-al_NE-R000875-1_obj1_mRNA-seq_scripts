#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce Euler diagrams of the treatment differentially expressed gene
# (DEG) list comparisons.
#-------------------------------------------------------------------------------
# Inputs:
# Results of the treatment DEG list comparisons.

# Outputs:
# .svg file of the Euler diagrams.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(eulerr)  # euler()
library(ggplotify)  # as.ggplot()
library(ggplot2)
library(ggpubr)  # ggarrange()

# LOAD DATA ----
# NOTE: This script assumes that the working directory is the 01_scripts 
# subdirectory.

# Change Directory ----

setwd("../02_outputs/00_all_tissues/")

treatmentComparisonResults <- 
  read.csv("22_NER0008751_obj1_sup_file_S1g_comparison_stats.csv")

# FUNCTION DEFINITIONS ----

MakeEulerPlot <- function (x, y) {
  # Make a Euler plot of the results of comparing the DEGs between control and 
  # removal treatments for a given tissue and direction of differential 
  # expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in time point 2 compared to time 
  #      point one) or down-regulated genes (genes more highly expressed in time
  #      point 1 compared to time point 2) ("up" or "down").
  #
  # Returns:
  #   ggplot2 object containing a Euler diagram of the results.
  
  # Assign variable based on arguments.
  
  if (x == "brain" && y == "up") {
    z <- 1
  } else if (x == "brain" && y == "down") {
    z <- 2
  } else if (x == "fatbody" && y == "up") {
    z <- 3
  } else if (x == "fatbody" && y == "down") {
    z <- 4
  } else if (x == "ovaries" && y == "up") {
    z <- 5
  } else if (x == "ovaries" && y == "down") {
    z <- 6
  } else {
    stop('Arguments x and/or y are incorrect.')
  }
  
  # Fit Euler Diagram to the data.
  
  eulerData <- euler(c(Removal = treatmentComparisonResults[z, 7],
                       Control = treatmentComparisonResults[z, 6],
                       "Removal&Control" = treatmentComparisonResults[z, 5]))
  
  # Plot Euler diagram.
  
  ePlot <- as.ggplot(plot(eulerData, 
                          fills = c("orange", "gray68", "white"),
                          quantities = TRUE))
  
  # Return Euler diagram.
  
  return(ePlot)
  
} 

# EXECUTED STATEMENTS ----

# Generate Euler Diagram of Each Comparison ----

# Brain.

brainUpPlot <- MakeEulerPlot("brain", "up")

brainDownPlot <- MakeEulerPlot("brain", "down")

# Fat body.

fatBodyUpPlot <- MakeEulerPlot("fatbody", "up")

fatBodyDownPlot <- MakeEulerPlot("fatbody", "down")

# Ovaries.

ovariesUpPlot <- MakeEulerPlot("ovaries", "up")

ovariesDownPlot <- MakeEulerPlot("ovaries", "down")

# Combine Diagrams into Single Plot ----

plotsCombined <- ggarrange(brainUpPlot, brainDownPlot,
                           fatBodyUpPlot, fatBodyDownPlot,
                           ovariesUpPlot, ovariesDownPlot,
                           labels = c("a", "b", "c", "d", "e", "f"),
                           ncol = 2, nrow = 3)

# Saving Combined Plot ----

# Save figures as SVG files so that the positioning of the labels
# can be adjusted manually for clarity.

ggsave("24_NER0008751_obj1_fig_5_euler_diagram.svg",
       plot = plotsCombined, height = 29.7, width = 21, units = "cm")
