#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 06.05.2021
# File last updated: 12.01.2023
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Produce a custom heatmap.
#-------------------------------------------------------------------------------
# Inputs:
# A matrix of gene names and expression values.

# Outputs:
# A heatmap with annotations.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(pheatmap)  # pheatmap().
library(ggplot2)  # ggsave(). 

# FUNCTION DEFINITION ----

PlotCustomHeatmap <- function (x) {
  # Plots a custom heat map, depending on the script calling the function.
  #
  # Args:
  #   x: number denoting which script is calling the function (41, 44 or 46).
  #
  # Returns:
  #   A heatmap.
  
  # Set column gaps based on argument x.
  
  if (x == 41) {
    columnGaps <- c(3, 6)
  } else if (x == 44 || x == 46) {
    columnGaps <- c(2, 4)
  }
  
  # Set column and row annotations based on argument x.
  
  if (x == 41) {
    colAnnotations <- 
      data.frame(Treatment = c(NA, "Removal queens", "Control queens", NA, 
                               "Removal queens", "Control queens", 
                               "Removal queens", "Control queens"),
                 Tissue = c("Brain", "Brain", "Brain", "Fat body", "Fat body", 
                            "Fat body", "Ovaries", "Ovaries"),
                 Species = c("D. melanogaster", "B. terrestris", "B. terrestris",
                             "D. melanogaster", "B. terrestris", "B. terrestris",
                             "B. terrestris", "B. terrestris"))
    
    row.names(colAnnotations) <- colnames(plotMatrix)
    
    rowAnnotations <- 
      data.frame(Longevity_influence = noBterZeroResults$longevity.influence)
    
    rownames(rowAnnotations) <- noBterZeroResults$name
    
    annotationColours <- list(Treatment = c("Removal queens" = "orange", "Control queens" = "gray68"),
                              Tissue = c("Brain" = "pink1", "Fat body" = "tan3", "Ovaries" = "white"),
                              Species = c("D. melanogaster" = "firebrick1", "B. terrestris" = "yellow"),
                              Longevity_influence =c("Pro-Longevity" = "white", "Anti-Longevity" = "black"))
  } else if (x == 44) {
    colAnnotations <- 
      data.frame(Treatment = c("Removal queens", "Control queens", "Removal queens", 
                               "Control queens", "Removal queens", "Control queens"),
                 Tissue = c("Brain", "Brain", "Fat body", "Fat body", 
                            "Ovaries", "Ovaries"))
    
    row.names(colAnnotations) <- colnames(plotMatrix)
    
    annotationColours <- list(Treatment = c("Removal queens" = "orange", "Control queens" = "gray68"),
                              Tissue = c("Brain" = "pink1", "Fat body" = "tan3", 
                                         "Ovaries" = "white"))
  } else if (x == 46) {
    colAnnotations <- 
      data.frame(Treatment = c("Removal queens", "Control queens", 
                               "Removal queens", "Control queens", 
                               "Removal queens", "Control queens"),
                 Tissue = c("Brain", "Brain", "Fat body", "Fat body", 
                            "Ovaries", "Ovaries"))
    
    row.names(colAnnotations) <- colnames(plotMatrix)
    
    annotationColours <- list(Treatment = c("Removal queens" = "orange", "Control queens" = "gray68"),
                              Tissue = c("Brain" = "pink1", "Fat body" = "tan3", 
                                         "Ovaries" = "white"))
  } else {
    stop('Argument x must equal 41, 44 or 46.')
  }
  
  # Plot heatmap based on argument x.
  
  if (x == 41) {
    heatmapPlot <- pheatmap(plotMatrix,
                            scale = "none",
                            cluster_rows = TRUE,
                            cluster_cols = FALSE,
                            annotation_col = colAnnotations,
                            annotation_row = rowAnnotations,
                            annotation_colors = annotationColours,
                            legend = TRUE,
                            cellwidth=15,
                            cellheight=10,
                            fontsize = 7, 
                            gaps_col = columnGaps)
  } else if (x == 44 || x == 46) {
    heatmapPlot <- pheatmap(plotMatrix,
                            scale = "none",
                            cluster_rows = TRUE,
                            cluster_cols = FALSE,
                            annotation_col = colAnnotations,
                            annotation_colors = annotationColours,
                            legend = TRUE,
                            cellwidth=15,
                            cellheight=10,
                            fontsize = 7, 
                            gaps_col = columnGaps)
  }
  
  # Return heatmap.
  
  return(heatmapPlot)
  
}
