#!/bin/bash
# -----------------------------------------------------------------------------
# Script: OTUs_analysis.sh
# Purpose: Process demultiplexed sequences with cutadapt trimming, merging, 
#          quality filtering, dereplication, contaminant removal, and OTU clustering.
# -----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero exit status.
set -e

FORWARD_PRIMER="CCTACGGGNGGCWGCAG"
REVESRE_PRIMER="GACTACHVGGGTATCTAATCC"

# # -----------------------------------------------------------------------------
# # OTUs: Trim Paired-End Sequences using Cutadapt
# # -----------------------------------------------------------------------------
# echo "Step 1 of 9: OTUs - Trimming paired-end sequences with cutadapt..."
# # Trim the paired-end reads using cutadapt.
# # Uses 60 cores, removes forward and reverse primer sequences,
# # discards sequences that do not contain the primer, and outputs trimmed sequences.
# qiime cutadapt trim-paired \
#   --i-demultiplexed-sequences $DATAFLOW/01-import_and_qc/paired-end-demux.qza \
#   --p-cores 60 \
#   --p-front-f $FORWARD_PRIMER \
#   --p-front-r $REVESRE_PRIMER \
#   --p-discard-untrimmed \
#   --o-trimmed-sequences demux-trimmed.qza \
#   --verbose
# echo "Step 1 of 9: OTUs - Cutadapt trimming completed."

# # -----------------------------------------------------------------------------
# # OTUs: Merge Paired-End Reads using VSEARCH
# # -----------------------------------------------------------------------------
# echo "Step 2 of 9: OTUs - Merging paired-end reads with VSEARCH..."
# # Merge the paired-end reads from the trimmed sequences.
# # Sets a minimum overlap of 15, allows up to 7 differences,
# # truncates reads at quality 15, and uses 30 threads.
# qiime vsearch merge-pairs \
#   --i-demultiplexed-seqs demux-trimmed.qza \
#   --o-merged-sequences vsearch-merged-seqs.qza \
#   --o-unmerged-sequences vsearch-unmerged-seqs.qza \
#   --p-minovlen 15 \
#   --p-maxdiffs 10 \
#   --p-truncqual 15 \
#   --p-threads 60 \
#   --verbose
# echo "Step 2 of 9: OTUs - Merging of paired-end reads completed."

# # -----------------------------------------------------------------------------
# # OTUs: Quality Filtering on Merged Sequences
# # -----------------------------------------------------------------------------
# echo "Step 3 of 9: OTUs - Quality filtering merged sequences..."
# # Apply a quality filter to the merged sequences using q-score filtering.
# qiime quality-filter q-score \
#   --i-demux vsearch-merged-seqs.qza \
#   --o-filtered-sequences filtered-merged-seqs.qza \
#   --o-filter-stats filtered-merged-stats.qza
# echo "Step 3 of 9: OTUs - Quality filtering completed."

# # -----------------------------------------------------------------------------
# # OTUs: Dereplicate Sequences using VSEARCH
# # -----------------------------------------------------------------------------
# echo "Step 4 of 9: OTUs - Dereplicating filtered sequences..."
# # Dereplicate the filtered merged sequences to generate a feature table and sequences.
# qiime vsearch dereplicate-sequences \
#   --i-sequences filtered-merged-seqs.qza \
#   --o-dereplicated-table derep-table.qza \
#   --o-dereplicated-sequences derep-seqs.qza \
#   --verbose
# echo "Step 4 of 9: OTUs - Dereplication completed."

# # -----------------------------------------------------------------------------
# # OTUs: Contaminant Identification using Prevalence Method
# # -----------------------------------------------------------------------------
# echo "Step 5 of 9: OTUs - Identifying contaminants in the dereplicated table..."
# # Identify contaminants based on the prevalence method using metadata and negative
# # control (H2O) as control indicator.
# qiime quality-control decontam-identify  \
#   --i-table derep-table.qza \
#   --m-metadata-file $METADATA \
#   --p-prev-control-column "horse" \
#   --p-prev-control-indicator "H2O" \
#   --p-method prevalence \
#   --verbose \
#   --o-decontam-scores decontam-scores.qza
# echo "Step 5 of 9: OTUs - Contaminant identification completed."

# # -----------------------------------------------------------------------------
# # OTUs: Remove Contaminants
# # -----------------------------------------------------------------------------
# echo "Step 6 of 9: OTUs - Removing contaminants from the dereplicated data..."
# # Remove features identified as contaminants from the dereplicated sequences.
# qiime quality-control decontam-remove \
#   --i-decontam-scores decontam-scores.qza \
#   --i-table derep-table.qza \
#   --i-rep-seqs derep-seqs.qza \
#   --p-threshold 0.05 \
#   --o-filtered-table decontam_table.qza \
#   --o-filtered-rep-seqs decontam_rep-seqs.qza


# Filter the feature table and representative sequences to remove positive and negative control
qiime feature-table filter-samples \
  --i-table decontam_table.qza \
  --m-metadata-file "$METADATA" \
  --p-where "[sample-id] IN ('PK', 'NK')" \
  --p-exclude-ids \
  --o-filtered-table filtered-table.qza

# Filter the representative sequences to remove positive and negative control
qiime feature-table filter-seqs \
  --i-data decontam_rep-seqs.qza \
  --i-table filtered-table.qza  \
  --o-filtered-data filtered-rep-seqs.qza

echo "Step 6 of 9: OTUs - Contaminant and controls removal completed."


# -----------------------------------------------------------------------------
# OTUs: OTU Clustering (De Novo) using VSEARCH
# -----------------------------------------------------------------------------
echo "Step 7 of 9: OTUs - Clustering features into OTUs at 97% identity..."
# Cluster sequences into OTUs at 97% sequence identity.
qiime vsearch cluster-features-de-novo \
  --i-table filtered-table.qza\
  --i-sequences filtered-rep-seqs.qza \
  --p-perc-identity 0.97 \
  --p-threads 60 \
  --o-clustered-table table-dn-97.qza \
  --o-clustered-sequences rep-seqs-dn-97.qza
echo "Step 7 of 9: OTUs - OTU clustering completed."

# -----------------------------------------------------------------------------
# OTUs: Summarize OTU Table and Visualize Representative Sequences
# -----------------------------------------------------------------------------
echo "Step 8 of 9: OTUs - Generating summary visualization for OTU table..."
# Generate a summary visualization of the clustered OTU table with metadata.
qiime feature-table summarize \
  --i-table table-dn-97.qza \
  --m-sample-metadata-file $METADATA \
  --o-visualization table-dn-97.qzv
echo "Step 8 of 9: OTUs - OTU table summary visualization completed."

echo "Step 9 of 9: OTUs - Visualizing representative sequences for OTUs..."
# Generate a visualization of the representative sequences.
qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dn-97.qza \
  --o-visualization rep-seqs-dn-97.qzv
echo "Step 9 of 9: OTUs - Representative sequence visualization completed."

echo "OTUs analysis completed successfully."
