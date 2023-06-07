#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce Euler diagrams of the all samples vs. no virus vs. with virus 
# differentially expressed gene (DEG) list comparisons.
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

# Load Comparison Statistics ----

# Change working directory.

setwd("../02_outputs/00_all_tissues/")

# Load data.

virusComparisonResults <- 
  read.csv("25_NER0008751_obj1_all_samples_vs_with_without_virus_comparison_stats.csv")

# Format data.

virusComparisonResults$Number_of_list_1_DEGs <- as.numeric(virusComparisonResults$Number_of_list_1_DEGs)

virusComparisonResults$Number_of_list_2_DEGs <- as.numeric(virusComparisonResults$Number_of_list_2_DEGs)

virusComparisonResults$Number_overlapping_DEGs <- as.numeric(virusComparisonResults$Number_overlapping_DEGs)

# FUNCTION DEFINITIONS ----

CompareThreeGeneLists <- function (x, y = "R_TP2G vs. C_TP2G", z) {
  # Compares overlap between three lists of differentially expressed genes 
  # (DEGs) to generate data for a Euler diagram.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain"
  #      "ovaries").
  #   y: string denoting the comparison between time points to be compared 
  #      (default = "R_TP2G vs. C_TP2G".)
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in the earlier chronological time 
  #      point compared to the later chronological time point) or down-regulated
  #      genes (genes more highly expressed in the later chronological time point
  #      compared to the earlier chronological time point) ("up" or "down").
  #
  # Returns:
  #   A data.frame of the results of the comparison.
  # NOTE: The working directory needs to be set to 02_outputs/0x_tissue/33_DESeq2_additional_DEG_list 
  #       corresponding to the tissue used for argument x.
  
  # Set variable based on argument x.
  
  if (x == "brain") {
    supFileLetter <- "r"
  } else if (x == "fatbody") {
    supFileLetter <- "s"
  } else if (x == "ovaries") {
    supFileLetter <- "t"
  } else {
    stop ('Argument x must equal "brain", "fatbody" or "ovaries".')
  }
  
  # Load data based on argument x.
  
  allSamplesDEGs <- read.csv(paste0("00_NER0008751_obj1_sup_file_S1", 
                                    supFileLetter, "_",
                                    x, "_all_samples_additional_DEG_list.csv"))
  
  noVirusDEGs <- read.csv(paste0("01_NER0008751_obj1_", x,
                                 "_no_virus_DEG_list.csv"))
  
  withVirusDEGs <- read.csv(paste0("02_NER0008751_obj1_", x,
                                   "_with_virus_DEG_list.csv"))
  
  # Index the data based on arguments y and z.
  
  compAllSamplesDEGs <- allSamplesDEGs[allSamplesDEGs[, "comparison"] == y, ]
  
  specificAllSamplesDEGs <- 
    compAllSamplesDEGs[compAllSamplesDEGs[, "expression_with_respect_to_chrono_age"] == paste0(z, "-regulated"), ]
  
  compNoVirusDEGs <- noVirusDEGs[noVirusDEGs[, "comparison"] == y, ]
  
  specificNoVirusDEGs <- 
    compNoVirusDEGs[compNoVirusDEGs[, "expression_with_respect_to_chrono_age"] == paste0(z, "-regulated"), ]
  
  compWithVirusDEGs <- withVirusDEGs[withVirusDEGs[, "comparison"] == y, ]
  
  specificWithVirusDEGs <- 
    compWithVirusDEGs[compWithVirusDEGs[, "expression_with_respect_to_chrono_age"] == paste0(z, "-regulated"), ]
  
  # Compare the 3 gene lists in pairs.
  
  allSamplesANDNoVirus <- intersect(specificAllSamplesDEGs$symbol, specificNoVirusDEGs$symbol)
  
  allSamplesANDWithVirus <- intersect(specificAllSamplesDEGs$symbol, specificWithVirusDEGs$symbol)
  
  noVirusANDWithVirus <- intersect(specificNoVirusDEGs$symbol, specificWithVirusDEGs$symbol)
  
  # Compare which genes are shared by all methods.
  
  sharedByAllThree <- intersect(allSamplesANDNoVirus, allSamplesANDWithVirus)
  
  # Determine genes shared by 2 out of 3 methods.
  
  allSamplesANDNoVirusNOTWithVirus <- setdiff(allSamplesANDNoVirus, sharedByAllThree)
  
  allSamplesANDWithVirusNOTNoVirus <- setdiff(allSamplesANDWithVirus, sharedByAllThree)
  
  noVirusANDWithVirusNOTAllSamples	<- setdiff(noVirusANDWithVirus, sharedByAllThree)
  
  # Determine genes only present in 1 of the methods.
  
  allSamplesOnly <- length(specificAllSamplesDEGs$symbol) - length(allSamplesANDNoVirusNOTWithVirus) - length(allSamplesANDWithVirusNOTNoVirus) - length(sharedByAllThree)
  
  noVirusOnly <- length(specificNoVirusDEGs$symbol) - length(allSamplesANDNoVirusNOTWithVirus) - length(noVirusANDWithVirusNOTAllSamples) - length(sharedByAllThree)
  
  withVirusOnly <- length(specificWithVirusDEGs$symbol) - length(allSamplesANDWithVirusNOTNoVirus) - length(noVirusANDWithVirusNOTAllSamples) - length(sharedByAllThree)
  
  # Collect data for euler diagram and output to file.
  
  eulerData <- 
    data.frame("category" = c("all_samples_only", "no_virus_only",
                              "with_virus_only", "all_samples&no_virus",
                              "all_samples&with_virus",
                              "no_virus&with_virus",
                              "all_samples&no_virus&with_virus"),
               "gene" = c(allSamplesOnly, noVirusOnly, withVirusOnly,
                          length(allSamplesANDNoVirusNOTWithVirus),
                          length(allSamplesANDWithVirusNOTNoVirus),
                          length(noVirusANDWithVirusNOTAllSamples),
                          length(sharedByAllThree)))
  
  # Return eulerData.
  
  return(eulerData)
  
}

MakeTwoWayEulerPlot <- function (x, y, z) {
  # Make a Euler plot of the results of comparing the DEGs between all samples 
  # and no virus analysis for a given tissue, treatment/time point and direction
  # of differential expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting the comparison between time points to be compared 
  #      ("R_TP1G vs. R_TP2G" or "R_TP2G vs. C_TP2G".)
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in the earlier chronological time 
  #      point compared to the later chronological time point) or down-regulated
  #      genes (genes more highly expressed in the later chronological time point
  #      compared to the earlier chronological time point) ("up" or "down").
  #
  # Returns:
  #   ggplot2 object containing a Euler diagram of the results.
  
  # Index data based on arguments.
  
  tissueData <- virusComparisonResults[virusComparisonResults[, "Tissue"] == x, ]
  
  comparisonData <- tissueData[tissueData[, "Comparison"] == y, ]
  
  a <- which(comparisonData$Direction_of_expression == paste0(z, "-regulated"))
  
  # Fit Euler Diagram to the data.
  
  eulerData <- euler(c("All Samples" = comparisonData[a, "Number_of_list_1_DEGs"] - comparisonData[a, "Number_overlapping_DEGs"],
                       "No virus" = comparisonData[a, "Number_of_list_2_DEGs"] - comparisonData[a, "Number_overlapping_DEGs"],
                       "All Samples&No virus" = comparisonData[a, "Number_overlapping_DEGs"]))
  
  # Plot Euler diagram.
  
  ePlot <- as.ggplot(plot(eulerData, 
                          fills = c("gray86", "white", "gray95"),
                          quantities = TRUE))
  
  # Return Euler diagram.
  
  return(ePlot)
  
} 

MakeThreeWayEulerPlot <- function (x, y = "R_TP2G vs. C_TP2G", z) {
  # Make a Euler plot of the results of comparing the DEGs between all samples, 
  # no virus samples and with virus samples for a given tissue, treatment/time 
  # point and direction of differential expression.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain"
  #      or "ovaries").
  #   y: string denoting the comparison between time points to be compared 
  #      ("R_TP2G vs. C_TP2G".)
  #   z: string denoting whether the DEG lists to compare are for up-regulated 
  #      genes (genes more highly expressed in the earlier chronological time 
  #      point compared to the later chronological time point) or down-regulated
  #      genes (genes more highly expressed in the later chronological time point
  #      compared to the earlier chronological time point) ("Up" or "Down").
  #
  # Returns:
  #   ggplot2 object containing a Euler diagram of the results.
  
  # Select data.
  
  tissueData <- eval(parse(text = paste0(x, z, "Comparison")))

  # Fit Euler Diagram to the data.
  
  eulerData <- euler(c("All Samples" = tissueData[1, "gene"],
                       "No virus" = tissueData[2, "gene"],
                       "With virus" = tissueData[3, "gene"],
                       "All Samples&No virus" = tissueData[4, "gene"],
                       "All Samples&With virus" = tissueData[5, "gene"],
                       "No virus&With virus" = tissueData[6, "gene"],
                       "All Samples&No virus&With virus" = tissueData[7, "gene"]))
  
  # Plot Euler diagram.
  
  ePlot <- as.ggplot(plot(eulerData, 
                          fills = c("gray86", "white", "gray28", "gray95", "gray45",
                                    "gray29", "gray60"),
                          quantities = TRUE))
  
  # Return Euler diagram.
  
  return(ePlot)
  
} 

# EXECUTED STATEMENTS ----

# Compare Gene Lists for Brain and Ovaries R_TP2G vs. C_TP2G ----

# Set working directory.

setwd("../02_brain/33_DESeq2_additional_DEG_list/")

# Compare lists.

brainUpComparison <- CompareThreeGeneLists("brain", z = "up")

brainDownComparison <- CompareThreeGeneLists("brain", z = "down")

# Reset working directory.

setwd("../../01_ovaries/33_DESeq2_additional_DEG_list/")

# Compare lists.

ovariesUpComparison <- CompareThreeGeneLists("ovaries", z = "up")

ovariesDownComparison <- CompareThreeGeneLists("ovaries", z = "down")

# Generate Euler Diagram of Each Comparison ----

# Brain.

brainRTP1vsRTP2UpPlot <- MakeTwoWayEulerPlot("brain", "R_TP1G vs. R_TP2G", "up")

brainRTP1vsRTP2DownPlot <- MakeTwoWayEulerPlot("brain", "R_TP1G vs. R_TP2G", "down")

brainRTP2vsCTP2UpPlot <- MakeThreeWayEulerPlot("brain", z = "Up")

brainRTP2vsCTP2DownPlot <- MakeThreeWayEulerPlot("brain", z = "Down")

# Fat body.

fatBodyRTP1vsRTP2UpPlot <- MakeTwoWayEulerPlot("fatbody", "R_TP1G vs. R_TP2G", "up")

fatBodyRTP1vsRTP2DownPlot <- MakeTwoWayEulerPlot("fatbody", "R_TP1G vs. R_TP2G", "down")

fatBodyRTP2vsCTP2UpPlot <- MakeTwoWayEulerPlot("fatbody", "R_TP2G vs. C_TP2G", "up")

fatBodyRTP2vsCTP2DownPlot <- MakeTwoWayEulerPlot("fatbody", "R_TP2G vs. C_TP2G", "down")

# Ovaries.

ovariesRTP1vsRTP2UpPlot <- MakeTwoWayEulerPlot("ovaries", "R_TP1G vs. R_TP2G", "up")

ovariesRTP1vsRTP2DownPlot <- MakeTwoWayEulerPlot("ovaries", "R_TP1G vs. R_TP2G", "down")

ovariesRTP2vsCTP2UpPlot <- MakeThreeWayEulerPlot("ovaries", z = "Up")

ovariesRTP2vsCTP2DownPlot <- MakeThreeWayEulerPlot("ovaries", z = "Down")

# Combine Diagrams into Single Plot ----

plotsCombined <- ggarrange(brainRTP1vsRTP2UpPlot, brainRTP1vsRTP2DownPlot,
                           brainRTP2vsCTP2UpPlot, brainRTP2vsCTP2DownPlot,
                           fatBodyRTP1vsRTP2UpPlot, fatBodyRTP1vsRTP2DownPlot,
                           fatBodyRTP2vsCTP2UpPlot, fatBodyRTP2vsCTP2DownPlot,
                           ovariesRTP1vsRTP2UpPlot, ovariesRTP1vsRTP2DownPlot,
                           ovariesRTP2vsCTP2UpPlot, ovariesRTP2vsCTP2DownPlot,
                           labels = c("a", "b", "c", "d", "e", "f", "g", "h",
                                      "i", "j", "k", "l"),
                           ncol = 2, nrow = 6)

# Saving Combined Plot ----

# Set working directory.

setwd("../../00_all_tissues/")

# Save figures as SVG files so that the positioning of the labels
# can be adjusted manually for clarity.

ggsave("27_NER0008751_obj1_fig_S19_all_samples_vs_with_without_virus_euler_diagram.svg",
       plot = plotsCombined, height = 27.7, width = 19, units = "cm")
