#!/bin/bash
#
#-------------------------------------------------------------------------------
# Author: David Prince
# Project: NER0008751 (Obj1) 
# Analysis: mRNA-seq
# Subsection: Pre-analysis tasks.
# Tasks: Download required Bombus terrestris (Bter), Apis mellifera and 
# Drosophila melanogaster files from Ensembl, and HoloBee file repositories.
#-------------------------------------------------------------------------------
# Inputs:
# None
#
# Outputs:
# Ensembl files and unzipped Holobee_v2016.1 archive.
#-------------------------------------------------------------------------------

# STEP 1: DOWNLOAD FASTA FILES ----
# NOTE: script assumes the working directory is 01_scripts.

# Change Directory ----

cd ../00_data/01_fasta

# Download Bombus terrestris Genome Fasta File and Checksums ----

# Download genome fasta file.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/Bombus_terrestris.Bter_1.0.dna.toplevel.fa.gz

# Download associated CHECKSUMS.
# Args:
# -O: Name of output file.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/dna/CHECKSUMS -O ../04_sums/Bter_dna_toplevel_checksums.txt

# Download Bombus terrestris Transcriptome Fasta File and Checksums ----

# Download transcriptome fasta file.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/cdna/Bombus_terrestris.Bter_1.0.cdna.all.fa.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/fasta/bombus_terrestris/cdna/CHECKSUMS -O ../04_sums/Bter_cdna_all_checksums.txt

# Download Bombus terrestris Peptide Sequences and Checksums ----

# Download peptide fasta file.

wget ftp://ftp.ensemblgenomes.org/pub/release-49/metazoa/fasta/bombus_terrestris/pep/Bombus_terrestris.Bter_1.0.pep.all.fa.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-49/metazoa/fasta/bombus_terrestris/pep/CHECKSUMS -O ../04_sums/Bter_pep_checksums.txt

# Download Apis mellifera Peptide Sequences and Checksums ----

# Download peptide fasta file.

wget ftp://ftp.ensemblgenomes.org/pub/release-49/metazoa/fasta/apis_mellifera/pep/Apis_mellifera.Amel_HAv3.1.pep.all.fa.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-49/metazoa/fasta/apis_mellifera/pep/CHECKSUMS -O ../04_sums/Amel_pep_checksums.txt

# Download Drosophila melanogaster Peptide Sequences and Checksums ----

# Download peptide fasta file.

wget ftp://ftp.ensemblgenomes.org/pub/release-49/metazoa/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.28.pep.all.fa.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-49/metazoa/fasta/drosophila_melanogaster/pep/CHECKSUMS -O ../04_sums/Dmel_pep_checksums.txt

# STEP 2: DOWNLOAD BTER ANNOTATION FILES ----

# Download GTF File and Checksums ----

# Change directory.

cd ../02_gtf

# Download GTF file.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/Bombus_terrestris.Bter_1.0.47.gtf.gz

# Download associated CHECKSUMS.

wget ftp://ftp.ensemblgenomes.org/pub/release-47/metazoa/gtf/bombus_terrestris/CHECKSUMS -O ../04_sums/Bter_gtf_checksums.txt

# STEP 3: DOWNLOAD BTER FEATURE TABLE ----

# Change Directory ----

cd ../05_feature_table

# Download Feature Table File ----

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/30195/102/GCF_000214255.1_Bter_1.0/GCF_000214255.1_Bter_1.0_feature_table.txt.gz

# Unzip features table.
# Args:
# --stdout: Keep file unchanged and write output to standard output.

gunzip --stdout GCF_000214255.1_Bter_1.0_feature_table.txt.gz > Bter_v1_feature_table.txt

# Download Associated CHECKSUMS ----

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/30195/102/GCF_000214255.1_Bter_1.0/md5checksums.txt -O ../04_sums/md5sums_Bter_feature_table.txt

# STEP 4: DOWNLOAD AND UNZIP HOLOBEE FILES ----

# Change Working Directory ----

cd ../01_fasta

# Download Holobee Archive ----

wget https://data.nal.usda.gov/system/files/HB_v2016.1.zip

# Note:
# Checksums for the Holobee data were downloaded separately from the website
# above and included in the file "holobee_md5.txt".

# Unzip Holobee Archive ----
# Args:
# --stdout: Keep file unchanged and write output to standard output.

unzip HB_v2016.1.zip

# Move Holobee Barcode File ----

mv HB_v2016.1/HB_Bar_v2016.1.fasta ./

# Remove Files no Longer Required ----
# Args:
# -r: remove directories and their contents recursively
# -f: do not ask permission to delete the files

rm -rf HB_v2016.1.zip

rm -rf __MACOSX

rm -rf HB_v2016.1
