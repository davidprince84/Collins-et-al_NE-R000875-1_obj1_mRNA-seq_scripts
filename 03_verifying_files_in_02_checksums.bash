#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Select relevant entries from checksums files and compare with the
# checksums generated from the downloaded files.
#-------------------------------------------------------------------------------
# Inputs:
# Insect fasta and annotation files, and associated checksum files, downloaded
# from Ensembl using script 02. Holobee fasta and checksums files downloaded
# using script 02.
#
# Outputs:
# .txt files stating whether the checksums match between the downloaded files
# and the original files.
#-------------------------------------------------------------------------------

# STEP 1: SELECT RELEVANT CHECKSUMS ----

# NOTE: This script assumes that the working directory is the 01_scripts subdirectory.

# Print Relevant Checksums to New Document for Comparison ----

# Change directory.

cd ../00_data/04_sums

# Print the relevant checksums.

# Bombus terrestris genome checksums.

awk '/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz/' Bter_dna_toplevel_checksums.txt > original_ensembl_checksums.txt

# B. terrestris cDNA checksums. 

awk '/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz/' Bter_cdna_all_checksums.txt >> original_ensembl_checksums.txt

# B. terrestris GTF checksums.

awk '/Bombus_terrestris.Bter_1.0.47.gtf.gz/' Bter_gtf_checksums.txt >> original_ensembl_checksums.txt
 
# B. terrestris peptide checksums.
 
awk '/Bombus_terrestris.Bter_1.0.pep.all.fa.gz/' Bter_pep_checksums.txt >> original_ensembl_checksums.txt

# A. mellifera peptide checksums.
 
awk '/Apis_mellifera.Amel_HAv3.1.pep.all.fa.gz/' Amel_pep_checksums.txt >> original_ensembl_checksums.txt
 
# D melanogaster peptide checksums.
 
awk '/Drosophila_melanogaster.BDGP6.28.pep.all.fa.gz/' Dmel_pep_checksums.txt >> original_ensembl_checksums.txt

# STEP 2: CONDUCT CHECKSUMS ON DOWNLOADED FILES ----
# Print the checksums results and the file name to a new document so that they can be compared.

# Conduct Checksums ----

# B. terrestris genome checksums.

echo "$(sum ../01_fasta/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz) $(echo Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz)" > downloaded_ensembl_checksums.txt

# B. terrestris cDNA checksums.

echo "$(sum ../01_fasta/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz) $(echo Bombus_terrestris.Bter_1.0.cdna.all.fa.gz)" >> downloaded_ensembl_checksums.txt

# B. terrestris GTF checksums.

echo "$(sum ../02_gtf/Bombus_terrestris.Bter_1.0.47.gtf.gz) $(echo Bombus_terrestris.Bter_1.0.47.gtf.gz)" >> downloaded_ensembl_checksums.txt

# B. terrestris peptide checksums.
 
echo "$(sum ../01_fasta/Bombus_terrestris.Bter_1.0.pep.all.fa.gz) $(echo Bombus_terrestris.Bter_1.0.pep.all.fa.gz)" >> downloaded_ensembl_checksums.txt
 
# A. mellifera peptide checksums.
  
echo "$(sum ../01_fasta/Apis_mellifera.Amel_HAv3.1.pep.all.fa.gz) $(echo Apis_mellifera.Amel_HAv3.1.pep.all.fa.gz)" >> downloaded_ensembl_checksums.txt
 
# D. melanogaster peptide checksums.
 
echo "$(sum ../01_fasta/Drosophila_melanogaster.BDGP6.28.pep.all.fa.gz) $(echo Drosophila_melanogaster.BDGP6.28.pep.all.fa.gz)" >> downloaded_ensembl_checksums.txt

# STEP 3: COMPARE THE ORIGINAL AND DOWNLOADED ENSEMBL CHECKSUMS ----
# Check that the checksums of the downloaded file match those provided with the file
# by comparing the files.

# Compare Checksums ----

# Make empty file for results output.

touch ../../02_outputs/00_all_tissues/00_checksums_ensembl_comparison_results.txt

# Compare results.

# grep arguments:
# -q: quiet, exits with zero status (mapped to true in the "if" statement) if any matches
# found, and exits with non-zero status (i.e. "false") if no matches found.
# -f fileName: obtains pattern(s) for comparison from fileName. 

while IFS= read -r FILE
do
if grep -q "$FILE" original_ensembl_checksums.txt
then
	RESULTS="$FILE \t Checksums match, therefore downloaded file is likely to be complete"
else 
	RESULTS="$FILE \t Checksums do not match, therefore downloaded file is likely to be incomplete"
fi
echo -e $RESULTS >> ../../02_outputs/00_all_tissues/00_checksums_ensembl_comparison_results.txt
done < "downloaded_ensembl_checksums.txt"

# STEP 4: COMPARE HOLOBEE MD5SUMS ----

# Print the Relevant MD5sums ----

awk '/HB_Bar_v2016.1.fasta/' md5sums_holobee.txt > original_md5sums_holobee.txt

# Copy and Compare the MD5sums ----

# Move a copy of the Holobee md5sums to the fasta directory.

cp original_md5sums_holobee.txt ../01_fasta

# Change directory.

cd ../01_fasta

# Compare the Holobee md5sums with the downloaded files and produce a text file of the results.

md5sum -c original_md5sums_holobee.txt > ../../02_outputs/00_all_tissues/01_md5sums_comparison_results.txt

# Remove copy of md5sums from fastq directory.

rm -f original_md5sums_holobee.txt

# STEP 5: COMPARE FEATURE TABLE MD5SUMS ----

# Print and Format the Relevant MD5sums ----

# Change directory.

cd ../04_sums

# Print md5sums.

awk '/GCF_000214255.1_Bter_1.0_feature_table.txt.gz/' md5sums_Bter_feature_table.txt > original_md5sums_feature_table.txt 

# Format md5sums by removing two characters before the file name (./).
# Pattern for sed is any two characters followed by G, replace with G.

sed 's/..G/G/' original_md5sums_feature_table.txt > formatted_md5sums_feature_table.txt

# Copy and Compare the MD5sums ----

# Move a copy of the md5sums to the fasta directory.

cp formatted_md5sums_feature_table.txt ../05_feature_table

# Change to the fasta directory.

cd ../05_feature_table

# Compare the md5sums with the downloaded files and produce a text file of the results.

md5sum -c formatted_md5sums_feature_table.txt >> ../../02_outputs/00_all_tissues/01_md5sums_comparison_results.txt

# Remove copy of md5sums from fastq directory.

rm -f formatted_md5sums_feature_table.txt
