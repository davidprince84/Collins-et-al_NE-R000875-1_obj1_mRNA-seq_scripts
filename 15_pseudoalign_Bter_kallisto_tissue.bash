#!/bin/bash

#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Quality control (QC) of the raw mRNA-seq reads.
# Tasks: Pseudoalign mRNA-seq reads to Bombus terrestris (Bter) transcriptome 
# using Kallisto. 
#-------------------------------------------------------------------------------
# Inputs:
# Raw mRNA-seq files. Bter transcriptome index for Kallisto.
#
# Outputs:
# Directories containing results of pseudoalignment (abundances).
# MultiQC summary of pseudoalignments.
# Text file containing study design table.
#-------------------------------------------------------------------------------

# LOAD SOFTWARE MODULES ----
# This is necessary when performing the analysis on a high performance computing
# cluster (HPC).  

module add kallisto/0.46.1
module add python/anaconda/2019.10/3.7  # Where MultiQC is installed.

# STEP 1: SET VARIABLES BASED ON ARGUMENT ----

# $1 takes the value of the first argument after the bash script name when the script is executed.

if [ ${1} = ovaries ]
then
	NUMBER=01
	TISSUE=ovary
elif [ ${1} = brain ]
then 
	NUMBER=02
	TISSUE=brain
elif [ ${1} = fatbody ]
then 
	NUMBER=03
	TISSUE=fatbody
else
	echo "Arguments incorrect, must be ovaries, brain or fatbody"
fi

# STEP 2: MAKE DIRECTORIES FOR KALLISTO OUTPUTs ----

mkdir ../02_outputs/${NUMBER}_${1}/22_kallisto_pseudoalignment_abundances

mkdir ../02_outputs/${NUMBER}_${1}/23_kallisto_pseudoalignment_summaries

# STEP 3: PSEUDOALIGN MRNA-SEQ READS USING KALLISTO ----

# Change Directory ----

cd ../00_data/00_fastq

# Pseudoalign with Kallisto ----

# Looping over mRNA-seq files for kallisto pseudoalignment.
# Args:
# -i: fileName for index to be used for quanitification.
# -o: directory to write output to.
# --rf-stranded: strand specific reads, first read reverse.
# NOTE: 2> captures standard error to file, which is where the pseudoalignment stats are sent.

for FILE in *${TISSUE}*R1.fastq.gz 
do
if [ ${1} = brain ]
then
	OUTPUT=${FILE:25:26}
	FILEHANDLE=${FILE:0:51}
elif [ ${1} = fatbody ]
then 
	OUTPUT=${FILE:25:28}
	FILEHANDLE=${FILE:0:53}
else [ ${1} = ovaries ]
	OUTPUT=${FILE:25:5}_${1}_${FILE:37:14}
	FILEHANDLE=${FILE:0:51}
fi
kallisto quant -i ../01_fasta/Bter_v1_cDNA_index.idx -o ../../02_outputs/${NUMBER}_${1}/22_kallisto_pseudoalignment_abundances/${OUTPUT} --rf-stranded ${FILE} ${FILEHANDLE}_R2.fastq.gz 2> ../../02_outputs/${NUMBER}_${1}/23_kallisto_pseudoalignment_summaries/${OUTPUT}.txt
done

# STEP 4: SUMMARISE PSEUDOALIGNMENTS WITH MULTIQC ----

# Change Directory ----

cd ../../02_outputs/${NUMBER}_${1}/23_kallisto_pseudoalignment_summaries/

# Run MultiQC ----
# -n fileName: call report fileName rather than default of multiqc_report.html
# -o: save report in directory 

multiqc . -n 14_NER0008751_obj1_${1}_kallisto_multiqc_report.html -o ../02_quality_control_reports

# STEP 5: GENERATE STUDY DESIGN FILE ----

# Initiate the file.
# Args:
# -e: enable interpretation of backslash escapes.

echo -e "sample\ttreatment\ttime_point\tlane" > ${1}_study_design.txt

# Loop over file names to fill in the details.

for NAME in DC*.txt 
do
if [ ${1} = brain ]
then
	SAMPLE=${NAME:12:11}
	TREATMENT=${NAME:12:1}
	TIMEPOINT=${NAME:14:4}
	LANE=${NAME:24:2}
elif [ ${1} = fatbody ]
then 
	SAMPLE=${NAME:14:11}
	TREATMENT=${NAME:14:1}
	TIMEPOINT=${NAME:16:4}
	LANE=${NAME:26:2}
else 
	SAMPLE=${NAME:14:11}
	TREATMENT=${NAME:14:1}
	TIMEPOINT=${NAME:16:4}
	LANE=${NAME:26:2}
fi
echo -e "${SAMPLE}\t${TREATMENT}\t${TIMEPOINT}\t${LANE}" >> ${1}_study_design.txt
done
