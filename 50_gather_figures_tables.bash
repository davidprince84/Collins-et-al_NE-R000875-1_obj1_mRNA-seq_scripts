#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq and comparative analysis with other species.
# Tasks: Make a new directory and copy all figure and supplementary files across.
#-------------------------------------------------------------------------------
# Inputs:
# Completed analysis (scripts #01 - #47)
#
# Outputs:
# A directory containing all the figures and supplementary files from the analysis.
#-------------------------------------------------------------------------------

# STEP 1: MAKE THE NEW DIRECTORY ----
# NOTE: script assumes the working directory is 01_scripts.

# Make Directory ----

mkdir ../02_outputs/20_all_figures_tables

# STEP 2: COPY RELEVANT FILES TO NEW DIRECTORY ----

# All Tissues ----

# Change directory.

cd ../02_outputs/00_all_tissues

# Copy figures.

cp *obj1_fig_* ../20_all_figures_tables

# Copy supplementary tables.

cp *obj1_table_S* ../20_all_figures_tables

# Ovaries ----

# Change directory.

cd ../01_ovaries/02_quality_control_reports

# Copy additional file.

cp *obj1_additional* ../../20_all_figures_tables

# Change directory.

cd ../31_DESeq2_DEG_lists

# Copy supplementary table. 

cp *obj1_table_S* ../../20_all_figures_tables

# Change directory.

cd ../32_DEG_top50_heatmaps

# Copy figure.

cp *obj1_fig_S* ../../20_all_figures_tables

# Brain ----

# Change directory.

cd ../../02_brain/02_quality_control_reports

# Copy additional file.

cp *obj1_additional* ../../20_all_figures_tables

# Change directory.

cd ../31_DESeq2_DEG_lists

# Copy supplementary table. 

cp *obj1_table_S* ../../20_all_figures_tables

# Change directory.

cd ../32_DEG_top50_heatmaps

# Copy figure.

cp *obj1_fig_S* ../../20_all_figures_tables

# Fat body ----

# Change directory.

cd ../../03_fatbody/02_quality_control_reports

# Copy additional file.

cp *obj1_additional* ../../20_all_figures_tables

# Change directory.

cd ../31_DESeq2_DEG_lists

# Copy supplementary table. 

cp *obj1_table_S* ../../20_all_figures_tables

# Change directory.

cd ../32_DEG_top50_heatmaps

# Copy figure.

cp *obj1_fig_S* ../../20_all_figures_tables

# STEP 3: REMOVE NUMBERS PRE-FIXING EACH FILE ----

# Change Directory ----

cd ../../20_all_figures_tables

# Remove Pre-fix Numbers ----

for FILE in *
do
	NEWFILE=${FILE:3}
mv $FILE $NEWFILE
done
