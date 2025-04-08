#!/bin/bash

# -----------------------------------------------------------------------------
# Script: run_quantification.sh
# Purpose: Perform phylogenetic diversity analysis, alpha & beta diversity 
#          assessments, taxonomy classification, and differential abundance testing.
# -----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero exit status.
set -e

# -----------------------------------------------------------------------------
# Step 1 - Set Variables and make Necessary Directories
# -----------------------------------------------------------------------------
DIR_NAME=$1
REP_SEQS=$2
TABLE=$3
SAMPLING_DEPTH=$4
METADATA_COLUMN=$5
METADATA_COLUMN2=$6
ALPHA_RARE_MAX_DEPTH=$7
SAMPLE_TYPE=$8
METADATA_COLUMN3=$9
REF_LEVEL=${10}
ANALYSIS_LEVEL=${11}
TAXA_COLLAPSE_LEVEL=${12}
TAG=${13}

echo "Step 1 of 12: Setting directories for analysis..."
mkdir -p $DATAFLOW/$DIR_NAME/{01-phylogeny,02-alpha_beta_diversity,03-taxonomy,04-differential_abundance,05-visualization}
cd "$DATAFLOW/$DIR_NAME/"
echo "Step 1 of 12: Directories set successfully."

# -----------------------------------------------------------------------------
# Step 2 - Generate Phylogenetic Tree
# -----------------------------------------------------------------------------
echo "Step 2 of 12: Generating a phylogenetic tree..."

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences $REP_SEQS \
  --p-n-threads 60 \
  --o-alignment 01-phylogeny/${TAG}_aligned-rep-seqs.qza \
  --o-masked-alignment 01-phylogeny/${TAG}_masked-aligned-rep-seqs.qza \
  --o-tree 01-phylogeny/${TAG}_unrooted-tree.qza \
  --o-rooted-tree 01-phylogeny/${TAG}_rooted-tree.qza
echo "Step 2 of 12: Phylogenetic tree generated."

export DIVERSITY=02-alpha_beta_diversity/diversity-core-metrics-phylogenetic

# -----------------------------------------------------------------------------
# Step 3 - Perform Core Diversity Analysis
# -----------------------------------------------------------------------------
echo "Step 3 of 12: Running core diversity analysis..."

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny 01-phylogeny/${TAG}_rooted-tree.qza \
  --i-table $TABLE \
  --p-n-jobs-or-threads 60 \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file $METADATA \
  --output-dir $DIVERSITY
echo "Step 3 of 12: Core diversity analysis completed."

# -----------------------------------------------------------------------------
# Step 4 - Perform Alpha Diversity Significance Tests
# -----------------------------------------------------------------------------
echo "Step 4 of 12: Calculating alpha diversity significance..."

qiime diversity alpha-group-significance \
  --i-alpha-diversity $DIVERSITY/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization $DIVERSITY/../${TAG}_faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity $DIVERSITY/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization $DIVERSITY/../${TAG}_evenness-group-significance.qzv

echo "Step 4 of 12: Alpha diversity significance analysis completed."

# -----------------------------------------------------------------------------
# Step 5 - Perform Beta Diversity Significance Tests
# -----------------------------------------------------------------------------
echo "Step 5 of 12: Calculating beta diversity significance..."

qiime diversity beta-group-significance \
  --i-distance-matrix $DIVERSITY/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column $METADATA_COLUMN \
  --p-pairwise \
  --o-visualization $DIVERSITY/../${TAG}_unweighted-unifrac-${METADATA_COLUMN}-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix $DIVERSITY/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column $METADATA_COLUMN2 \
  --p-pairwise \
  --o-visualization $DIVERSITY/../${TAG}_unweighted-unifrac-${METADATA_COLUMN2}-group-significance.qzv

echo "Step 5 of 12: Beta diversity significance analysis completed."

# -----------------------------------------------------------------------------
# Step 6 - Generate Alpha Rarefaction Curve
# -----------------------------------------------------------------------------
echo "Step 6 of 12: Generating alpha rarefaction curve..."

qiime diversity alpha-rarefaction \
  --i-table $TABLE \
  --i-phylogeny 01-phylogeny/${TAG}_rooted-tree.qza \
  --p-max-depth $ALPHA_RARE_MAX_DEPTH \
  --m-metadata-file $METADATA \
  --o-visualization $DIVERSITY/../${TAG}_alpha-rarefaction.qzv
echo "Step 6 of 12: Alpha rarefaction curve completed."

# -----------------------------------------------------------------------------
# Step 7 - Download and Apply Taxonomic Classifier
# -----------------------------------------------------------------------------
echo "Step 7 of 12: Downloading taxonomic classifier..."

CLASSIFIER="$CLASS_DIR/16SrRNA_V3-4_classifier.qza"


echo "Step 7 of 12: Taxonomic classifier downloaded."

echo "Step 8 of 12: Performing taxonomic classification..."

qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads $REP_SEQS \
  --o-classification 03-taxonomy/${TAG}_taxonomy.qza

qiime metadata tabulate \
  --m-input-file 03-taxonomy/${TAG}_taxonomy.qza \
  --o-visualization 03-taxonomy/${TAG}_taxonomy.qzv

qiime taxa barplot \
  --i-table $TABLE \
  --i-taxonomy 03-taxonomy/${TAG}_taxonomy.qza \
  --m-metadata-file $METADATA \
  --o-visualization 03-taxonomy/${TAG}_taxa-bar-plots.qzv

echo "Step 8 of 12: Taxonomic classification completed."

# -----------------------------------------------------------------------------
# Step 9 - Filter Samples for Differential Abundance Analysis
# -----------------------------------------------------------------------------
echo "Step 9 of 12: Filtering samples for differential abundance analysis..."

qiime feature-table filter-samples \
  --i-table $TABLE \
  --m-metadata-file $METADATA \
  --p-where "[$METADATA_COLUMN3] IN ('$REF_LEVEL','$ANALYSIS_LEVEL')" \
  --o-filtered-table 04-differential_abundance/${TAG}_sample_type-table.qza

echo "Step 9 of 12: Sample filtering completed."

# -----------------------------------------------------------------------------
# Step 10 - Run ANCOM-BC for Differential Abundance Analysis
# -----------------------------------------------------------------------------
echo "Step 10 of 12: Performing differential abundance analysis using ANCOM-BC..."

qiime composition ancombc \
  --i-table 04-differential_abundance/${TAG}_sample_type-table.qza \
  --m-metadata-file $METADATA \
  --p-formula $METADATA_COLUMN3 \
  --p-reference-levels $METADATA_COLUMN3::$REF_LEVEL \
  --o-differentials 04-differential_abundance/${TAG}_ancombc-condition.qza

qiime composition da-barplot \
  --i-data 04-differential_abundance/${TAG}_ancombc-condition.qza \
  --p-significance-threshold 0.001 \
  --o-visualization 04-differential_abundance/${TAG}_da-barplot-condition.qzv

echo "Step 10 of 12: Differential abundance analysis completed."

# -----------------------------------------------------------------------------
# Step 11 - Collapse Taxonomy and Perform ANCOM-BC on Collapsed Data
# -----------------------------------------------------------------------------
echo "Step 11 of 12: Collapsing taxonomy and performing ANCOM-BC analysis..."

qiime taxa collapse \
  --i-table 04-differential_abundance/${TAG}_sample_type-table.qza \
  --i-taxonomy 03-taxonomy/${TAG}_taxonomy.qza \
  --p-level $TAXA_COLLAPSE_LEVEL \
  --o-collapsed-table 04-differential_abundance/${TAG}_sample_type-table-l${TAXA_COLLAPSE_LEVEL}.qza

qiime tools export \
  --input-path 04-differential_abundance/${TAG}_sample_type-table-l${TAXA_COLLAPSE_LEVEL}.qza \
  --output-path 04-differential_abundance/${TAG}_export

biom convert -i 04-differential_abundance/${TAG}_export/feature-table.biom -o 05-visualization/${TAG}_genus.tsv --to-tsv

qiime composition ancombc \
  --i-table 04-differential_abundance/${TAG}_sample_type-table-l${TAXA_COLLAPSE_LEVEL}.qza \
  --m-metadata-file $METADATA \
  --p-formula $METADATA_COLUMN3 \
  --p-reference-levels $METADATA_COLUMN3::$REF_LEVEL \
  --o-differentials 04-differential_abundance/${TAG}_l${TAXA_COLLAPSE_LEVEL}-ancombc-condition.qza

qiime composition da-barplot \
  --i-data 04-differential_abundance/${TAG}_l${TAXA_COLLAPSE_LEVEL}-ancombc-condition.qza \
  --p-significance-threshold 0.001 \
  --o-visualization 04-differential_abundance/${TAG}_l${TAXA_COLLAPSE_LEVEL}-da-barplot-condition.qzv

echo "Step 11 of 12: Taxonomic collapse and differential abundance analysis completed."


# -----------------------------------------------------------------------------
# Step 12 - Run aldex2 plug-in for generating differentials.qza file needed in R visualisation
# -----------------------------------------------------------------------------

# NOTE: The output of this plug-in is needed for R visualisation (qiime2R) to
# generate e.g. the volcano plot and phylogenetic tree. Since permissions are not
# given to install this plug-in on the server, the scirpt will scip this step.

# echo "Step 12 of 12: - Running aldex2 plug-in for differential abundance analysis..."  


# qiime aldex2 aldex2 \
#     --i-table 04-differential_abundance/${TAG}_gum-table.qza \
#     --m-metadata-file $METADATA \
#     --m-metadata-column disease_state  \
#     --output-dir 05-visualization/gum-test \
#     --verbose

# echo "Step 12 of 12: - Running aldex2 plug-in for differential abundance analysis..."  


# echo "Step 12 of 12: - aldex2 plug-in analysis completed."


echo "Analysis pipeline completed successfully!"