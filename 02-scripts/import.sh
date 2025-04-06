#!/bin/bash
# -----------------------------------------------------------------------------
# Script: import.sh
# Purpose: Generate a manifest file for QIIME2 paired-end sequencing import,
#          import the fastq files into QIIME2, and generate quality visualizations.
# -----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero exit status.
set -e

# -----------------------------------------------------------------------------
# Step 1: Write the Manifest File
# -----------------------------------------------------------------------------

# This step is required to map sample identifiers to fast.gz or fastq absolute filepaths.
# The manifest file also indicates the direction of the reads in each file. 
# For more information how to make the manifest file visit teh qiime2 dcoumentation: https://docs.qiime2.org/2024.10/tutorials/importing/
echo "Writing manifest file for QIIME2 import..."

# Write the header for the manifest file (tab-delimited) as required by QIIME2.
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "$MANIFEST"

# Loop through each forward fastq file in the DATA directory.
for f1 in "$DATA"/*_1.fastq.gz; do
  # Extract the sample name by removing the trailing '_1.fastq.gz'
  sample_name=$(basename "${f1%_1.fastq.gz}")
  
  # Define the corresponding reverse fastq file by replacing '_1.fastq.gz' with '_2.fastq.gz'
  f2="$DATA/${sample_name}_2.fastq.gz"
  
  # Check if the reverse file exists; if yes, append the sample details to the manifest.
  if [[ -f "$f2" ]]; then
    echo -e "$sample_name\t$DATA/$(basename "$f1")\t$DATA/$(basename "$f2")" >> "$MANIFEST"
  else
    echo "Warning: Missing reverse file for sample '$sample_name'" >&2
  fi
done

echo "Manifest file stored at: $MANIFEST"

# -----------------------------------------------------------------------------
# Step 2: Import Paired-End Fastq Files into QIIME2
# -----------------------------------------------------------------------------
echo "Starting QIIME2 import of paired-end fastq files..."

# Import the paired-end sequencing data using the manifest file. In this case study
# we use demuliplexed paired-end fastq files including a positive (E.coli) and negative (water) control.
# The qiime2 import artificat will have the format "PairedEndFastqManifestPhred64V2"
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$MANIFEST" \
  --output-path paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

# Display metadata information about the imported artifact.
qiime tools peek paired-end-demux.qza

# Generate a summary visualization to assess the quality of the imported sequences.
qiime demux summarize \
    --i-data paired-end-demux.qza \
    --o-visualization paired-end_quality.qzv

echo "QIIME2 import was successful."
