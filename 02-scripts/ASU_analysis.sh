#!/bin/bash
# -----------------------------------------------------------------------------
# Script: ASU_analysis.sh
# Purpose: Process paired-end sequence data using QIIME2’s DADA2 plugin, 
#          generate quality control visualizations, and remove contaminants.
# -----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero exit status.
set -e

# -----------------------------------------------------------------------------
# Step 1 - DADA2 Denoising: Process paired-end reads to remove errors and generate a 
# feature table and representative sequences.
# -----------------------------------------------------------------------------
echo "Step 1 of 9: ASV - Starting DADA2 denoising..."

# Remove the first N bases (length of forward and reverse primers) and low
# quality regions based on the visualization of the quality scores in paired-end_quality.qzv:
# Truncate forward reads to 267 bases and reverse reads to 221 bases.
# Use 60 threads for parallel processing.
# Output feature table, representative sequences, and denoising statistics.
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs $DATAFLOW/01-import_and_qc/paired-end-demux.qza \
  --p-trim-left-f $LENGTH_F_PRIMER \
  --p-trim-left-r $LENGTH_R_PRIMER \
  --p-trunc-len-f 265 \
  --p-trunc-len-r 225 \
  --p-n-threads 60 \
  --o-table dada2-table.qza \
  --o-representative-sequences dada2-seqs.qza \
  --o-denoising-stats dada2-stats.qza \
  --verbose

echo "Step 1 of 9: ASV - DADA2 denoising completed."

# -----------------------------------------------------------------------------
# Step 2 - Visualize Denoising Statistics
# -----------------------------------------------------------------------------
echo "Step 2 of 9: Visualizing denoising statistics..."

# Generate a visualization of the DADA2 statistics.
qiime metadata tabulate \
    --m-input-file dada2-stats.qza \
    --o-visualization dada2-stats.qzv

echo "Step 2 of 9: ASV - Denoising statistics visualization completed."

# -----------------------------------------------------------------------------
# Step 3 - Summarize Feature Table
# -----------------------------------------------------------------------------
echo "Step 3 of 9: ASV - Summarizing feature table..."

# Generate a summary visualization of the feature table, including sample metadata.
qiime feature-table summarize \
    --i-table dada2-table.qza \
    --m-sample-metadata-file $METADATA \
    --o-visualization dada2-table.qzv 

echo "Step 3 of 9: ASV - Feature table summary completed."

# -----------------------------------------------------------------------------
# Step 4 - Visualize Representative Sequences
# -----------------------------------------------------------------------------
echo "Step 4 of 9: ASV - Visualizing representative sequences..."

# Make a visualization table of the representative sequences.
qiime feature-table tabulate-seqs \
    --i-data dada2-seqs.qza \
    --o-visualization dada2-seqs.qzv

echo "Step 4 of 9: ASV - Representative sequences visualization completed."

# -----------------------------------------------------------------------------
# Step 5 - Contaminant Identification
# -----------------------------------------------------------------------------
echo "Step 5 of 9: ASV - Identifying contaminants using prevalence method..."

# Identify potential contaminants using the prevalence method.
# Compares feature prevalence in negative control (H2O).
qiime quality-control decontam-identify \
  --i-table dada2-table.qza \
  --m-metadata-file $METADATA \
  --p-method prevalence \
  --p-prev-control-column "horse" \
  --p-prev-control-indicator "H2O" \
  --o-decontam-scores prev_decontam_scores.qza \
  --verbose

echo "Step 5 of 9: ASV - Contaminant identification completed."

# -----------------------------------------------------------------------------
# Step 6 - Visualize Contaminant Scores
# -----------------------------------------------------------------------------
echo "Step 6 of 9: ASV - Visualizing contaminant scores..."

# Generate a visualization of the decontamination scores 
# to assess which features may be contaminants.
qiime quality-control decontam-score-viz \
  --i-decontam-scores prev_decontam_scores.qza \
  --i-table dada2-table.qza \
  --i-rep-seqs dada2-seqs.qza \
  --p-threshold 0.1 \
  --p-no-weighted \
  --p-bin-size 0.05 \
  --o-visualization decontan_score_viz.qzv

echo "Step 6 of 9: ASV - Contaminant score visualization completed."

# -----------------------------------------------------------------------------
# Step 7 - Remove Contaminants
# -----------------------------------------------------------------------------
echo "Step 7 of 9: ASV - Removing contaminants from data..."

# Filter out features that exceed the contaminant threshold and removing rows of controls.
qiime quality-control decontam-remove \
  --i-decontam-scores prev_decontam_scores.qza \
  --i-table dada2-table.qza \
  --i-rep-seqs dada2-seqs.qza \
  --p-threshold 0.05 \
  --o-filtered-table decontam-table.qza \
  --o-filtered-rep-seqs decontam-rep-seqs.qza

# Filter the feature table and representative sequences to remove positive and negative control
qiime feature-table filter-samples \
  --i-table decontam-table.qza \
  --m-metadata-file "$METADATA" \
  --p-where "[sample-id] IN ('PK', 'NK')" \
  --p-exclude-ids \
  --o-filtered-table filtered-table.qza

# Filter the representative sequences to remove positive and negative control
qiime feature-table filter-seqs \
  --i-data decontam-rep-seqs.qza \
  --m-metadata-file "$METADATA" \
  --p-where "[sample-id] IN ('PK', 'NK')" \
  --p-exclude-ids \
  --o-filtered-data filtered-rep-seqs.qza

echo "Step 7 of 9: ASV - Contaminant removal completed."

# -----------------------------------------------------------------------------
# Step 8 - Summarize Filtered Feature Table
# -----------------------------------------------------------------------------
echo "Step 8 of 9: ASV - Summarizing filtered feature table..."

# Generate a summary visualization of the filtered feature table 
# using the provided sample metadata.
qiime feature-table summarize \
    --i-table filtered-table.qza \
    --m-sample-metadata-file $METADATA \
    --o-visualization filtered-table.qzv 

echo "Step 8 of 9: ASV - Filtered feature table summary completed."

# -----------------------------------------------------------------------------
# Step 9 - Visualize Filtered Representative Sequences
# -----------------------------------------------------------------------------
echo "Step 9 of 9: ASV - Visualizing filtered representative sequences..."

# Make a visualization table for the filtered representative sequences.
qiime feature-table tabulate-seqs \
    --i-data filtered-rep-seqs.qza \
    --o-visualization filtered-rep-seqs.qzv

echo "Step 9 of 9: ASV - Filtered representative sequences visualization completed."
echo "Analysis completed successfully!"
