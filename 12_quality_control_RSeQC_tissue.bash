#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Run RSeQC geneBody_coverage.py and junction_saturation.py. 
# Combine reports using MultiQC.
#-------------------------------------------------------------------------------
# Inputs:
# Sorted bam files of reads aligned with HISAT2. Indexes for sorted bam files.
#
# Outputs:
# MultiQC html report and associated files.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module add python/anaconda/2019.10/3.7  # Where RSeQC and MultiQC are installed.
module add R/4.0.0  # Needed for RSeQC for graphing.

# STEP 1: SET VARIABLES BASED ON ARGUMENTS ----

# $1 takes the value of the first argument after the bash script name when the script is executed.

if [ ${1} = ovaries ]
then
	NUMBER=01
elif [ ${1} = brain ]
then 
	NUMBER=02

elif [ ${1} = fatbody ]
then 
	NUMBER=03
else
	echo "Arguments incorrect, must be ovaries, brain or fatbody"
fi

# STEP 2: RSEQC ANALYSIS ----

# Change Directory ----

cd ../02_outputs/${NUMBER}_${1}/10_hisat2_aligned_reads/

# Make Results Directory ----

mkdir ../12_rseqc_results

# Run geneBody_coverage.py ----

# Visualizes coverage of mapped reads over length of genes.

# geneBody_coverage.py arguments:
# -r fileName: gene annotation file in BED12 format called fileName.
# -i directory/Path: a path to the input file directory. 
# -o fileName: prefix of the output files.

geneBody_coverage.py -r ../../../00_data/03_bed/bter_v1.BED -i ./ -o ../12_rseqc_results/RSeQC_${1}_results

# Run junction_saturation.py ----

# Visualizes the number of junctions discovered for a given read depth.

# Loop passes each file to the junction_saturation.py script and selects the output name.

# junction_saturation.py arguments.
# -i fileName: input file in bam format.
# -r file/Path/fileName: path and file name of gene annotation file in BED12 format.
# -o fileName: prefix for output files.

for FILE in *.bam
do
if [ ${1} = brain ]
then
	OUTFILE=${FILE:12:14}
elif [ ${1} = fatbody ]
then 
	OUTFILE=${FILE:14:14}
else [ ${1} = ovaries]
	OUTFILE=${FILE:12:14}
fi
junction_saturation.py -i ${FILE} -r ../../../00_data/03_bed/bter_v1.BED -o ../12_rseqc_results/${OUTFILE}
done

# STEP 3: MULTIQC OF RSEQC REPORTS ----

# Change Directory ----

cd ../12_rseqc_results

# Run MultiQC ----

# Run MultiQC python module to collate all alignment summary reports in folder together into one report
# -n fileName: call report fileName rather than default of multiqc_report.html.
# -o file/path: create report in subdirectory file/path.

multiqc . -n 12_NER0008751_obj1_${1}_rseqc_multiqc.html -o ../02_quality_control_reports
