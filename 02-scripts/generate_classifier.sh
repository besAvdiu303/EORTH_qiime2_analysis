#!/bin/bash

# -----------------------------------------------------------------------------
# Script: generate_classifier.sh
# Purpose: A classifier for the V3-V4 region of the 16S rRNA gene using the SILVA database
#          is needed since no classifier is available for this region in the QIIME2
# -----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero exit status.
set -e

echo "Generating a classifier for the V3-V4 region of the 16S rRNA gene using the SILVA database..."

# -----------------------------------------------------------------------------
# Step 1 - Set Variables and make Necessary Directories
# -----------------------------------------------------------------------------

# Set variables
F_PRIMER=$1
R_PRIMER=$2

# -----------------------------------------------------------------------------
# Step 2 - Download SILVA Database
# -----------------------------------------------------------------------------
# Download SILVA database version 138.2 (SSURef_NR99) including species labels
# Outputs:
#   silva_ref_seq.qza - Reference sequences
#   silva_ref_tax.qza - Reference taxonomy
echo "Step 1 started ..."
qiime rescript get-silva-data \
  --p-version "138.2" \
  --p-target "SSURef_NR99" \
  --p-include-species-labels \
  --o-silva-sequences silva_ref_seq.qza \
  --o-silva-taxonomy silva_ref_tax.qza \
  --verbose

# -----------------------------------------------------------------------------
# Step 3 - Cull Sequences
# -----------------------------------------------------------------------------
# Remove sequences that are duplicates or contain ambiguous bases
# Input: silva_ref_seq.qza
# Output: silva_ref_seq.qza (overwrites input with cleaned sequences)

echo "Step 2 started ..."
qiime rescript cull-seqs \
  --i-sequences silva_ref_seq.qza \
  --o-clean-sequences silva_ref_seq.qza \
  --p-n-jobs 60 \
  --verbose

# -----------------------------------------------------------------------------
# Step 4 - Filter Eukaryotic Sequences
# -----------------------------------------------------------------------------
# Remove eukaryotic sequences from the database
# Inputs:
#   silva_ref_seq.qza - Reference sequences
#   silva_ref_tax.qza - Reference taxonomy
# Output: silva_ref_seq.qza (overwrites input with filtered sequences)

echo "Step 3 started ..."
qiime taxa filter-seqs \
  --i-sequences silva_ref_seq.qza \
  --i-taxonomy silva_ref_tax.qza \
  --p-exclude "d__Eukaryota" \
  --p-mode "contains" \
  --o-filtered-sequences silva_ref_seq.qza \
  --verbose

# -----------------------------------------------------------------------------
# Step 5 - Filter Sequences by Length
# -----------------------------------------------------------------------------
# Filter sequences based on taxonomic classification and minimum length requirements
# Archaea: minimum 900 bp
# Bacteria: minimum 1200 bp
# Inputs:
#   silva_ref_seq.qza - Reference sequences
#   silva_ref_tax.qza - Reference taxonomy
# Outputs:
#   silva_ref_seq.qza - Filtered sequences meeting length requirements
#   silva_discard.qza - Discarded sequences

echo "Step 4 started ..."
qiime rescript filter-seqs-length-by-taxon \
  --i-sequences silva_ref_seq.qza \
  --i-taxonomy silva_ref_tax.qza \
  --p-labels "Archaea" "Bacteria" \
  --p-min-lens 900 1200 \
  --o-filtered-seqs silva_ref_seq.qza \
  --o-discarded-seqs silva_discard.qza \
  --verbose

# -----------------------------------------------------------------------------
# Step 6 - Extract Primer-Specific Reads
# -----------------------------------------------------------------------------
# Extract reads matching the provided primer sequences (V3-V4 region)
# Input: silva_ref_seq.qza
# Output: silva_extract_reads.qza - Extracted reads matching primers

echo "Step 5 started ..."
qiime feature-classifier extract-reads \
  --i-sequences silva_ref_seq.qza \
  --p-f-primer "$F_PRIMER" \
  --p-r-primer "$R_PRIMER" \
  --p-n-jobs 60 \
  --o-reads silva_extract_reads.qza \
  --verbose

# -----------------------------------------------------------------------------
# Step 7 - Dereplicate Sequences
# -----------------------------------------------------------------------------
# Remove duplicate sequences and corresponding taxonomy
# Inputs:
#   silva_extract_reads.qza - Extracted reads
#   silva_ref_tax.qza - Reference taxonomy
# Outputs:
#   silva_extract_reads.qza - Dereplicated sequences
#   silva_derep_taxa.qza - Dereplicated taxonomy

echo "Step 6 started ..."
qiime rescript dereplicate \
  --i-sequences silva_extract_reads.qza \
  --i-taxa silva_ref_tax.qza \
  --o-dereplicated-sequences silva_extract_reads.qza \
  --o-dereplicated-taxa silva_derep_taxa.qza \
  --p-threads 60 \
  --verbose

# -----------------------------------------------------------------------------
# Step 8 - Train Naive Bayes Classifier
# -----------------------------------------------------------------------------
# Train a Naive Bayes classifier using the dereplicated sequences and taxonomy
# Inputs:
#   silva_extract_reads.qza - Dereplicated sequences
#   silva_derep_taxa.qza - Dereplicated taxonomy
# Output: 16SrRNA_V3-4_classifier.qza - Trained classifier

echo "Step 7 started ..."
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads silva_extract_reads.qza \
  --i-reference-taxonomy silva_derep_taxa.qza \
  --o-classifier 16SrRNA_V3-4_classifier.qza \
  --verbose


echo "The 16SrRNA_V3-4_classifier is generated succesfully."