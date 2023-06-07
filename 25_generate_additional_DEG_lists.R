#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Differential expression analysis.
# Tasks: Produce additional lists of differentially expressed genes (DEGs) for 
# each tissue comparing between the treatments.
#-------------------------------------------------------------------------------
# Inputs:
# Bter_v1_transcripts2genes.txt file, genesymbols2genenames.txt file, Kallisto 
# abundances, virus_samples.csv

# Outputs:
# .csv files containing lists of DEGs for a given tissue and treatment.
# .csv file containing all genes expressed in the analysis, to be used as a
# background list in downstream comparisons.
#-------------------------------------------------------------------------------

# LOADING PACKAGES ----

# DESeq2 and tximport loaded via the scripts in the 
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

# Gene Symbol to Gene Name File ----

setwd("../05_feature_table/")

symbolToName <- read.csv("genesymbols2genenames.txt")

# Virus Samples ----

setwd("../11_virus_samples")

virus_samples <- read.csv(file = "virus_samples.csv",
                          col.names = c("sample", "virus"))

# LOAD CUSTOM FUNCTIONS ----

setwd("../../01_scripts")

source("f_BuildDDSForDEAnalysis.R")

# FUNCTION DEFINITIONS ----

listDEGsBetweenTreatments <- function (x, y, z, a = "all_samples") {
  # Generate a list of DEGs for a given tissue for two specified treatment/
  # time-point combinations.
  #
  # Args:
  #   x: string denoting which tissue is to be analysed ("brain", "fatbody" 
  #      or "ovaries").
  #   y: string denoting the name of first treatment/time-point (chronologically)
  #      to be compared ("R_TP1G", "R_TP2G" or "C_TP1G").
  #   z: string denoting the name of second treatment/time-point (chronologically)
  #      to be compared ("R_TP2G", "C_TP1G" or "C_TP2G").
  #   a: string denoting which set of samples should be analysed ("with_virus", 
  #      "no_virus" or "all_samples" (default)).
  #
  # Returns:
  #   Data.frame of DEGs.
  
  # Assign variables based on arguments.
  
  if (x == "brain") {
    ddsName <- "ddsBrain"
  } else if (x == "fatbody") {
    ddsName <- "ddsFatBody"
  } else if (x == "ovaries") {
    ddsName <- "ddsOvaries"
  } else {
    stop('Argument x must be "brain", "fatbody" or "ovaries".')
  }
  
  if (a == "all_samples") {
    samplesName <- "AllSamples"
  } else if (a == "no_virus") {
    samplesName <- "NoVirus"
  } else if (a == "with_virus") {
    samplesName <- "WithVirus"
  } else {
    stop('Argument a must be "with_virus", "no_virus" or "all_samples" (default).')
  }
  
  dds <- eval(parse(text = paste0(ddsName, samplesName)))
  
  # Extract differential expression results for chosen comparison.
  
  allExtractedResults <- 
    as.data.frame(results(dds, 
                          contrast = c("condition", z, y),
                          alpha = 0.05, 
                          lfcThreshold = 0))
  
  # Remove instances where padj = NA.
  
  allExtractedResults <- 
    allExtractedResults[!(is.na(allExtractedResults[, "padj"])), ]
  
  # Index the significantly differentially expressed genes into separate object.
  
  degResults <- 
    allExtractedResults[allExtractedResults[, "padj"] < 0.05, ]
  
  # Add column of gene symbols from the row names.
  
  degResults$symbol <- row.names(degResults)
  
  # Add gene names to results by merging with symbolToName.
  
  namedResults <- merge(degResults, symbolToName, all.x = TRUE)
  
  # Add column stating comparison.
  
  namedResults$comparison <- rep(paste0(y, " vs. ", z), 
                                 times = length(namedResults[, 1]))
  
  # Add column stating whether the gene is "up-regulated" (more expressed) or 
  # "down-regulated" (less expressed) with respect to chronological time.
  
  for (ROW in 1:length(namedResults[, 1])) {
    if (namedResults[ROW, "log2FoldChange"] > 0) {
      namedResults[ROW, "expression_with_respect_to_chrono_age"] <- "up-regulated"
    } else if (namedResults[ROW, "log2FoldChange"] < 0) {
      namedResults[ROW, "expression_with_respect_to_chrono_age"] <- "down-regulated"
    }
  }
  
  # Return namedResults.
  
  return(namedResults)
  
}

# EXECUTED STATEMENTS ----

# Build dds Objects for Each Tissue ----

ddsBrainAllSamples <- BuildDDSForDEAnalysis("brain")

ddsBrainNoVirus <- BuildDDSForDEAnalysis("brain", y = "no_virus")

ddsBrainWithVirus <- BuildDDSForDEAnalysis("brain", y = "with_virus")

ddsFatBodyAllSamples <- BuildDDSForDEAnalysis("fatbody")

ddsFatBodyNoVirus <- BuildDDSForDEAnalysis("fatbody", y = "no_virus")

ddsOvariesAllSamples <- BuildDDSForDEAnalysis("ovaries")

ddsOvariesNoVirus <- BuildDDSForDEAnalysis("ovaries", y = "no_virus")

ddsOvariesWithVirus <- BuildDDSForDEAnalysis("ovaries", y = "with_virus")

# Make Gene Lists ----

# Brain, all samples.

brainRTP1GvsCTP1GAllSamplesList <- listDEGsBetweenTreatments("brain", "R_TP1G", "C_TP1G")

brainRTP2GvsCTP1GAllSamplesList <- listDEGsBetweenTreatments("brain", "R_TP2G", "C_TP1G")

brainRTP2GvsCTP2GAllSamplesList <- listDEGsBetweenTreatments("brain", "R_TP2G", "C_TP2G")

brainRTP1GvsRTP2GAllSamplesList <- listDEGsBetweenTreatments("brain", "R_TP1G", "R_TP2G")

# Brain, no virus.

brainRTP1GvsRTP2GNoVirusList <- listDEGsBetweenTreatments("brain", "R_TP1G", "R_TP2G", a = "no_virus")

brainRTP2GvsCTP2GNoVirusList <- listDEGsBetweenTreatments("brain", "R_TP2G", "C_TP2G", a = "no_virus")

# Brain, with virus.

brainRTP2GvsCTP2GWithVirusList <- listDEGsBetweenTreatments("brain", "R_TP2G", "C_TP2G", a = "with_virus")

# Combine DEG lists.

combinedBrainAllSamplesLists <- rbind(brainRTP1GvsCTP1GAllSamplesList, brainRTP2GvsCTP1GAllSamplesList,
                                      brainRTP2GvsCTP2GAllSamplesList, brainRTP1GvsRTP2GAllSamplesList)

combinedBrainNoVirusLists <- rbind(brainRTP1GvsRTP2GNoVirusList,
                                   brainRTP2GvsCTP2GNoVirusList)

# Fat body, all samples.

fatBodyRTP1GvsCTP1GAllSamplesList <- listDEGsBetweenTreatments("fatbody", "R_TP1G", "C_TP1G")

fatBodyRTP2GvsCTP1GAllSamplesList <- listDEGsBetweenTreatments("fatbody", "R_TP2G", "C_TP1G")
# No differentially expressed genes returned.

fatBodyRTP2GvsCTP2GAllSamplesList <- listDEGsBetweenTreatments("fatbody", "R_TP2G", "C_TP2G")

fatBodyRTP1GvsRTP2GAllSamplesList <- listDEGsBetweenTreatments("fatbody", "R_TP1G", "R_TP2G")

# Fat body, no virus

fatBodyRTP1GvsRTP2GNoVirusList <- listDEGsBetweenTreatments("fatbody", "R_TP1G", "R_TP2G", a = "no_virus")

fatBodyRTP2GvsCTP2GNoVirusList <- listDEGsBetweenTreatments("fatbody", "R_TP2G", "C_TP2G", a = "no_virus")

# Combine DEG lists.

combinedFatBodyAllSamplesLists <- rbind(fatBodyRTP1GvsCTP1GAllSamplesList, #fatBodyRTP2GvsCTP1GList, 
                                        fatBodyRTP2GvsCTP2GAllSamplesList,
                                        fatBodyRTP1GvsRTP2GAllSamplesList)

combinedFatBodyNoVirusLists <- rbind(fatBodyRTP1GvsRTP2GNoVirusList,
                                     fatBodyRTP2GvsCTP2GNoVirusList)

# Ovaries, all samples.

ovariesRTP1GvsCTP1GAllSamplesList <- listDEGsBetweenTreatments("ovaries", "R_TP1G", "C_TP1G")

ovariesRTP2GvsCTP1GAllSamplesList <- listDEGsBetweenTreatments("ovaries", "R_TP2G", "C_TP1G")

ovariesRTP2GvsCTP2GAllSamplesList <- listDEGsBetweenTreatments("ovaries", "R_TP2G", "C_TP2G")

ovariesRTP1GvsRTP2GAllSamplesList <- listDEGsBetweenTreatments("ovaries", "R_TP1G", "R_TP2G")

# Ovaries, no virus.

ovariesRTP1GvsRTP2GNoVirusList <- listDEGsBetweenTreatments("ovaries", "R_TP1G", "R_TP2G", a = "no_virus")

ovariesRTP2GvsCTP2GNoVirusList <- listDEGsBetweenTreatments("ovaries", "R_TP2G", "C_TP2G", a = "no_virus")

# Ovaries, with virus.

ovariesRTP2GvsCTP2GWithVirusList <- listDEGsBetweenTreatments("ovaries", "R_TP2G", "C_TP2G", a = "with_virus")

# Combine DEG lists.

combinedOvariesAllSamplesLists <- rbind(ovariesRTP1GvsCTP1GAllSamplesList, 
                                      ovariesRTP2GvsCTP1GAllSamplesList,
                                      ovariesRTP2GvsCTP2GAllSamplesList,
                                      ovariesRTP1GvsRTP2GAllSamplesList)

combinedOvariesNoVirusLists <- rbind(ovariesRTP1GvsRTP2GNoVirusList,
                                   ovariesRTP2GvsCTP2GNoVirusList)

# Save Gene Lists ----

# Brain.

dir.create("../02_outputs/02_brain/33_DESeq2_additional_DEG_list")
# Will produce a warning if directory already exists.

# Change directory

setwd("../02_outputs/02_brain/33_DESeq2_additional_DEG_list")

# Save lists.

write.csv(combinedBrainAllSamplesLists, 
          "00_NER0008751_obj1_table_S18_brain_all_samples_additional_DEG_list.csv", 
          row.names = FALSE)

write.csv(combinedBrainNoVirusLists, 
          "01_NER0008751_obj1_brain_no_virus_DEG_list.csv", 
          row.names = FALSE)

write.csv(brainRTP2GvsCTP2GWithVirusList, 
          "02_NER0008751_obj1_brain_with_virus_DEG_list.csv", 
          row.names = FALSE)

# Fat body.

dir.create("../../03_fatbody/33_DESeq2_additional_DEG_list")
# Will produce a warning if directory already exists.

# Change directory

setwd("../../03_fatbody/33_DESeq2_additional_DEG_list")

# Save lists.

write.csv(combinedFatBodyAllSamplesLists, 
          "00_NER0008751_obj1_table_S19_fatbody_all_samples_additional_DEG_list.csv", 
          row.names = FALSE)

write.csv(combinedFatBodyNoVirusLists, 
          "01_NER0008751_obj1_fatbody_no_virus_DEG_list.csv", 
          row.names = FALSE)

# Ovaries.

dir.create("../../01_ovaries/33_DESeq2_additional_DEG_list")
# Will produce a warning if directory already exists.

# Change directory

setwd("../../01_ovaries/33_DESeq2_additional_DEG_list")

# Save lists.

write.csv(combinedOvariesAllSamplesLists, 
          "00_NER0008751_obj1_table_S20_ovaries_all_samples_additional_DEG_list.csv", 
          row.names = FALSE)

write.csv(combinedOvariesNoVirusLists, 
          "01_NER0008751_obj1_ovaries_no_virus_DEG_list.csv", 
          row.names = FALSE)

write.csv(ovariesRTP2GvsCTP2GWithVirusList, 
          "02_NER0008751_obj1_ovaries_with_virus_DEG_list.csv", 
          row.names = FALSE)
