#-------------------------------------------------------------------------------
# Author: David Prince
# File started: 22.03.2021
# File last updated: 13.01.2023
# Project: NER0008751 (Obj1) 
# Analysis: Comparative analysis with other species.
# Tasks: Compare Bombus terrestris (Bter) differentially 
# expressed genes (DEGs) with Drosophila melanogaster (Dmel) TI-J-LiFe network 
# genes, as defined by Korb et al. (2021), or Dmel enzymatic antioxidant network 
# genes, as defined by Kramer et al. (2021).
#-------------------------------------------------------------------------------
# Inputs:
# Bter DEGs from the current study. OrthoFinder results. List of Dmel
# TI-J-LiFe genes from Korb et al. (2021). List of Dmel enzymatic antioxidant 
# network genes from Kramer et al. (2021).

# Outputs:
# List with two objects: 
# 1) Data.frame recording all the results of the statistical tests.
# 2) Data.frame recording all the genes that overlap between the Dmel genes and 
#    the DEGs from the current study. 
#-------------------------------------------------------------------------------

# FUNCTION DEFINITION ----

CompareBterDEGsAndDmelGenes <- function (x, y, z = "all", a) {
  # Compares the overlap in single-copy orthologous genes between Bter DEG lists 
  # from the current study and Dmel genes in the TI-J-LiFe or enzymatic antioxidant 
  # network gene lists.
  #
  # Args:
  #   x: string denoting which tissue to compare DEG lists from ("brain",
  #      "fatbody" or "ovaries").
  #   y: string denoting which treatment the B. terrestris DEG list is from
  #      ("C" or "R").
  #   z: number denoting the number of DEGs from B. terrestris to compare with 
  #      the TI-J-LiFe or enzymatic antioxidant network genes. Number will 
  #      select the top z most positive and negative expressed DEGs based on 
  #      log fold change, therefore overall number of DEGs is 2 x z. 
  #      (Default = all DEGs).
  #   a: string denoting which Dmel gene list the Bter genes are being compared 
  #      with ("life" or "antioxidant").
  #
  # Returns:
  #   List with two objects: 
  #   1) Data.frame recording all the results of the statistical tests.
  #   2) Data.frame recording all the genes that overlap between the Dmel genes 
  #      and the DEGs from the current study. 
  # NOTE: The working directory needs to be set to the location of the gene 
  # list for the specified tissue.
  
  # Load Bter lists.
  
  bterUpDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                "up_regulated_LFC0.csv"))
  
  bterDownDEGs <- read.csv(paste0(x, "_results_", y, "_treatment_", 
                                  "down_regulated_LFC0.csv"))
  
  backgroundList <- read.csv(paste0(x, "_all_expressed_genes_and_DE_results.csv"))
  
  # Reduce number of columns and rename column in background list.
  
  backgroundList <- as.data.frame(backgroundList[, 1])
  
  colnames(backgroundList) <- "Bter_gene_ID"
  
  # Combine Bter up- and down-regulated genes.
  
  bterDEGs <- rbind(bterUpDEGs, bterDownDEGs)  
  
  # Filter DEGs, as specified by argument z.
  
  if (z == "all") {
    filteredBterDEGs <- bterDEGs  
  } else if ((2 * z) > length(bterDEGs[, 1])) {
    print('Argument z is too large, there are not suffient DEGs, using all DEGs instead.')
    filteredBterDEGs <- bterDEGs
    z <- "all, as z too large"
  } else {
    # Make all log2FoldChange values absolute, so that they can be ordered by
    # change irrespective of direction.
    
    bterDEGs$log2FoldChange <- abs(bterDEGs$log2FoldChange)
    
    # Sort data.frame so that rows are sorted by logFoldChange.
    
    sortedBterDEGs <- 
      bterDEGs[order(-(bterDEGs$log2FoldChange)), ]
    
    # Index the top 2 x z genes.
    
    filteredBterDEGs <- sortedBterDEGs[(1:(2 * z)), ]
    
  }
  
  # Simplify to just the gene symbols and rename column.
  
  filteredBterDEGs <- as.data.frame(filteredBterDEGs[, c("symbol", "name")])
  
  colnames(filteredBterDEGs) <- c("Bter_gene_ID", "Bter_gene_name")
  
  # Add Dmel orthologues.
  
  bterWithOrthos <- merge(filteredBterDEGs, dmelBterOrthologues)
  
  backgroundWithOrthos <- merge(backgroundList, dmelBterOrthologues)
  
  # Add orthologues to the Dmel list, based on argument a.
  
  if (a == "life") {
    dmelListOrthos <- merge(TIJLiFeGenes, dmelBterOrthologues)
    # 81 single-copy orthologues from the 123 Dmel genes.
    
    listID <- "TI-J-LiFe"
    
    maxLength <- 81
    
  } else if (a == "antioxidant") {
    dmelListOrthos <- merge(antioxidantGenes, dmelBterOrthologues)
    # 16 single-copy orthologues from 58 Dmel genes.
    
    listID <- "antioxidant genes"
    
    maxLength <- 16
    
  } else {
    stop('Argument a must equal "life" or "antioxidant.')
  }
  
  # Initialise results data.frame.
  
  resultsDF <- data.frame("Tissue" = x,
                          "Treatment" = y,
                          "Number_top_+/-_DEGs" = z,
                          "Dmel_list" = listID,
                          "Number_Bter_DEGs" = length(bterWithOrthos$Bter_gene_ID),
                          "Number_Dmel_genes" = length(dmelListOrthos$Flybase_gene_ID),
                          "Number_overlapping_genes" = NA,
                          "Number_genes_only_in_Bter" = NA,
                          "Number_genes_only_in_Dmel" = NA,
                          "Number_genes_in_neither" = NA,
                          "p_value" = NA,
                          "Odds_ratio" = NA,
                          "Alpha_value" = NA,
                          "Percentage_of_Dmel_genes_present" = NA) 
  
  # Calculate overlap between Bter and Dmel.
  
  # Determine overlapping genes and rename columns.
  
  overlappingGenes <- bterWithOrthos[bterWithOrthos$Bter_gene_ID %in% dmelListOrthos$Bter_gene_ID, ]
  
  valuesOfArgs <- paste(x, y, z, sep = "_")
  
  colnames(overlappingGenes) <- c(paste0(valuesOfArgs, "_Bter_gene_ID"),
                                  paste0(valuesOfArgs, "_Bter_gene_name"),
                                  "Orthogroup",
                                  paste0(valuesOfArgs, "_Flybase_gene_ID"))
   
  # Write results to data.frame, with extra blank space so that 
  # all results can be combined later.
  
  overlappingGenesDF <- overlappingGenes[, c(1, 2, 4)]
  
  if (is.na(overlappingGenesDF[1, 1])) {
    overlappingGenesDF[1, ] <- "No overlapping genes"
  }
  
  fillStart <- length(overlappingGenesDF[, 1]) + 1
  
  overlappingGenesDF[fillStart:maxLength, ] <- ""
  
  # Add results to data.frame.
  
  resultsDF$Number_overlapping_genes <- length(overlappingGenes[, 1])
  
  resultsDF$Percentage_of_TI_J_LiFe_present <- 
    round((resultsDF$Number_overlapping_genes/resultsDF$Number_Dmel_genes)*100, digits = 1)
  
  # Determine non-overlapping genes and add to data.frame.
  
  resultsDF$Number_genes_only_in_Bter <- 
    resultsDF$Number_Bter_DEGs - resultsDF$Number_overlapping_genes
  
  resultsDF$Number_genes_only_in_Dmel <-
    resultsDF$Number_Dmel_genes - resultsDF$Number_overlapping_genes
  
  # Determine number of non-DEGs and add to data.frame.
  
  resultsDF$Number_genes_in_neither <- 
    length(backgroundWithOrthos$Bter_gene_ID) - resultsDF$Number_overlapping_genes - 
    resultsDF$Number_genes_only_in_Bter - resultsDF$Number_genes_only_in_Dmel
  
  # Perform statistical test.
  
  # Create a matrix representing the numbers of genes in both lists, only 
  # Bter list, only Dmel list, and neither list.
  
  contingencyTable <- matrix(c(resultsDF$Number_overlapping_genes, 
                               resultsDF$Number_genes_only_in_Bter,
                               resultsDF$Number_genes_only_in_Dmel, 
                               resultsDF$Number_genes_in_neither))
  
  # Change the dimensions to 2 rows and 2 columns.
  
  dim(contingencyTable) <- c(2,2)
  
  # Conduct two-tailed Fisher's Exact Test on the results
  # to determine whether the number of shared genes between the two lists
  # is significantly higher or lower than expected by chance.
  
  fisherResults <- fisher.test(contingencyTable)
  
  # Add results of statistical tests to resultsDF.
  
  resultsDF$p_value <- fisherResults$p.value
  
  resultsDF$Odds_ratio <- fisherResults$estimate[[1]]
  
  # Return list of results.
  
  listToReturn <- list("genes" = overlappingGenesDF, 
                       "stats" = resultsDF)
  
  return(listToReturn)
  
}
