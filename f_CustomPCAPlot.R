#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 11.01.2021
# File last updated: 21.12.2022
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Tasks: Generate custom principal component analysis (PCA) plot.
#-------------------------------------------------------------------------------
# Inputs:
# dds (DESeq2 data set) object of Kallisto abundances.

# Outputs:
# ggplot2 object of the PCA plot.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

library(ggrepel)  # geom_text_repel()
# ggrepel loads ggplot2, needed for ggplot() and associated functions.

# FUNCTION DEFINITION ----

CustomPCAPlot <- function(x, y = "no", intgroup = "condition", ntop = 2000, techreps = "collapsed", virussymbols = "no") {
  # Acknowledgement:
  # This function was adapted from the pcaPlot function in the DESeq2 R package
  # v1.28.1.
  #
  # Makes a custom principal component analysis (PCA) plot from a dds object
  # and uses symbols to denote whether samples have many virus-aligning reads.
  #
  # Args:
  #   x: string denoting the tissue to produce a PCA plot for ("brain", 
  #      "fatbody" or "ovaries").
  #   y: string denoting whether samples containing many virus-aligning
  #      reads should be excluded from the analysis ("yes" and "no" (default)).
  #   intgroup: interesting groups: a character vector of names in colData(x) to 
  #             use for grouping (default of "condition").
  #   ntop: number of top genes to use for principal components, selected by 
  #         highest row variance (default of 2000)
  #   techreps: string denoting whether the technical replicates in the sequencing
  #             have been collapsed by DESeq2 ("collapsed" (default) and 
  #             "uncollapsed").
  #   virussymbols: string denoting whether samples containing many virus-aligning
  #                reads should be represented by a different symbols ("no" 
  #                (default) or "yes"). 
  #                Note: to use this property of the function
  #                the file "virus_samples.csv" must have been loaded into the
  #                environment.
  #
  # Returns:
  #   ggplot2 PCA plot.
  
  # Set variables based on argument.
  
  if (x == "brain") {
    ddsName <- "ddsBrain"
  } else if (x == "fatbody") {
    ddsName <- "ddsFatBody"
  } else if (x == "ovaries") {
    ddsName <- "ddsOvaries"
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
  
  # Transform the data to stabilize the variance across the mean.
  # Select  rlog transformation chosen if there are relatively few 
  # samples (n < 30), or vst transformation if there are more than 30
  # samples.
  
  if (length(dds$sample) < 30) {
    ddsTransform <- rlog(dds)  
  } else if (length(dds$sample) >= 30) {
    ddsTransform <- vst(dds)
  }
  
  # Calculate the variance for each gene.
  
  rv <- rowVars(assay(ddsTransform))
  
  # Select the ntop genes by variance.
  
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  
  # Perform a PCA on the data in assay(x) for the selected genes.
  
  pca <- prcomp(t(assay(ddsTransform)[select, ]))
  
  # The contribution to the total variance for each component.
  
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  if (!all(intgroup %in% names(colData(ddsTransform)))) {
    stop("The argument 'intgroup' should specify columns of colData(dds).")
  }
  
  intgroup.df <- as.data.frame(colData(ddsTransform)[, intgroup, drop=FALSE])
  
  # Add the intgroup factors together to create a new grouping factor.
  
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(ddsTransform)[[intgroup]]
  }
  
  # Determine whether name labels are for collapsed or uncollapsed technical
  # replicates.
  
  if (techreps == "collapsed") {
    nameLabels <- colnames(ddsTransform)
  } else if (techreps == "uncollapsed") {
    nameLabels <- paste0(ddsTransform$sample, "_", ddsTransform$lane)
  } else {
    stop('Argument techreps must be either "collapsed" (default) or "uncollapsed".')
  }
  
  # Assemble the data for the plot.
  
  d <- data.frame(PC1 = pca$x[,1], 
                  PC2 = pca$x[,2], 
                  group = group, 
                  intgroup.df, 
                  name = nameLabels)
  
  # Include values for virus-aligning samples if desired, and generate plot.
  
  if (virussymbols == "yes") {
    d2 <- merge(d, virus_samples, by.x = "name")
    pcaPlot <- ggplot(data = d2, aes_string(x = "PC1", y = "PC2", color = "group", shape = "virus"))
  } else if (virussymbols == "no") {
    pcaPlot <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group"))
  }
  
  # Add labels to plot.
  
  pcaPlotLabels <- pcaPlot +
    geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() +
    geom_text_repel(aes(label = name), size = 2, max.overlaps = Inf)
  
  # Format plot.
  
  formattedPCAPlot <- pcaPlotLabels + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Return plot.
  
  return(formattedPCAPlot)
  
}
