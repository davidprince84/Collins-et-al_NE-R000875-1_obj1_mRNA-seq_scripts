#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Perform FastQC analysis, combine the results in each tissue into a 
# report with MultiQC, and produce a summary text file of elements to check 
# further.
#-------------------------------------------------------------------------------
# Inputs:
# Tissue specific argument from command line. mRNA-seq files in fastq.gz format. 
#
# Outputs:
# FastQC results. MultiQC report. Text files containing the "fails" and 
# "warnings" from FastQC.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC). 

module add fastqc/0.11.9
module add python/anaconda/2019.10/3.7  # Where MultiQC module is installed.

# STEP 1: SET VARIABLES BASED ON ARGUMENT ----

# $1 takes the value of the first argument after the bash script name when the script is executed.

if [ ${1} = ovaries ]
then
	NUMBER=01
	SUPFIG=5
	TISSUE=ovary
elif [ ${1} = brain ]
then 
	NUMBER=02
	SUPFIG=3
	TISSUE=brain
elif [ ${1} = fatbody ]
then 
	NUMBER=03
	SUPFIG=4
	TISSUE=fatbody
else
	echo "Arguments incorrect, must be ovaries, brain or fatbody"
fi

# STEP 2: RUNNING FASTQC ----
# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Change Directory ----

cd ../00_data/00_fastq

# Run FastQC on all mRNA-seq Files ----

fastqc *${TISSUE}*.fastq*

# STEP 3: MAKING DIRECTORIES FOR RESULTS ----

# Making Directories ----

# Making results directories in tissue-specific output directories.

mkdir ../../02_outputs/${NUMBER}_${1}/00_fastqc_raw_reads

mkdir ../../02_outputs/${NUMBER}_${1}/01_fastqc_summaries_raw_reads

mkdir ../../02_outputs/${NUMBER}_${1}/02_quality_control_reports

# STEP 4: PARSING FASTQC RESULTS ----

# Moving Results ----

# Move all FastQC output files (.html and .zip) to appropriate tissue specific outputs folder.

mv *${TISSUE}*.zip ../../02_outputs/${NUMBER}_${1}/00_fastqc_raw_reads

mv *${TISSUE}*.html ../../02_outputs/${NUMBER}_${1}/00_fastqc_raw_reads

# Creating Reports of FastQC Data ----

# Change working directory to the relevant outputs folder.

cd ../../02_outputs/${NUMBER}_${1}/00_fastqc_raw_reads

# Run MultiQC python module to collate all FastQC reports in folder together into one report.
# -n fileName: call report fileName rather than default of multiqc_report.html.
# -o file/path: create report in subdirectory file/path.

multiqc . -n 00_NER0008751_obj1_additional_file_${SUPFIG}_${1}_fastqc_multiqc.html -o ../02_quality_control_reports

# Unzip all the FastQC results in folder using a for loop.

for FILENAME in *.zip
do
unzip $FILENAME
done

# Concatenate all the summary.txt files from the FastQC reports.

cat */summary.txt > ../01_fastqc_summaries_raw_reads/${1}_fastqc_combined_summaries.txt

# Change working directory to 10_FastQC_summaries.

cd ../01_fastqc_summaries_raw_reads

# Generate documents listing FastQC "FAIL" and "WARN" results.

grep FAIL ${1}_fastqc_combined_summaries.txt > ${1}_fastqc_fails.txt

grep WARN ${1}_fastqc_combined_summaries.txt > ${1}_fastqc_warnings.txt

# Sort reports by column 2 so that categories are grouped together and produce a new file.

sort -k 2 -o 01_${1}_fastqc_combined_summaries_sorted.txt ${1}_fastqc_combined_summaries.txt

sort -k 2 -o 02_${1}_fastqc_fails_sorted.txt ${1}_fastqc_fails.txt 

sort -k 2 -o 03_${1}_fastqc_warnings_sorted.txt ${1}_fastqc_warnings.txt
