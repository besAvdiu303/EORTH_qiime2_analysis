#!/bin/bash
# -----------------------------------------------------------------------------
# Script: run.sh
# Purpose: Automate data processing, ASV/OTUs analysis, and visualization with logging
# -----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero exit status.
set -e

# -----------------------------------------------------------------------------
# Initial Setup: Define directories, environment variables, and logging
# -----------------------------------------------------------------------------

echo "================================================================================"
echo "                    Starting QIIME2 Data Processing Pipeline                           "
echo "================================================================================"

# Set the working directory to the current directory. It must be the directory where this script is stored
export RUNPATH="$PWD"

# Define paths to key directories and files.
export DATA="$RUNPATH/01-data/20241209-raw_data"
export SCRIPTS="$RUNPATH/02-scripts"
export DATAFLOW="$RUNPATH/03-dataflow"
export MANIFEST="$DATAFLOW/01-import_and_qc/manifest.csv"

# Define the metadata file path. 
# The original metadata was validated using the open source Google sheets add-on "Keemei"
# For further information, refer to:
# Keemei: cloud-based validation of tabular bioinformatics file formats in Google Sheets. 
# Rideout JR, Chase JH, Bolyen E, Ackermann G, González A, Knight R, 
# Caporaso JG. GigaScience. 2016;5:27. http://dx.doi.org/10.1186/s13742-016-0133-6
export METADATA="$DATA/../metadata.tsv"
export CLASS_DIR="$DATA/../classifier/"
mkdir -p "$CLASS_DIR"

# Calculate and store the length of the primers for the V3-V4 region of the 16S rRNA gene.
F_PRIMER="CCTACGGGNGGCWGCAG"
R_PRIMER="GACTACHVGGGTATCTAATCC"

export LENGTH_F_PRIMER=$(echo -n $F_PRIMER | wc -c)
export LENGTH_R_PRIMER=$(echo -n $R_PRIMER | wc -c)

# -----------------------------------------------------------------------------
# Directory Structure Setup
# Make necessary directories for raw data and dataflow if they don't exist.
# -----------------------------------------------------------------------------

mkdir -p "$DATA"  # Ensure raw data directory exists.
mkdir -p "$RUNPATH/03-dataflow"/{01-import_and_qc,02-ASV_denosing,04-OTUs_clustering}

# Make a log directory under RUNPATH.
LOG_DIR="$RUNPATH/04-log"
mkdir -p "$LOG_DIR"

# Define the log file.
LOG_FILE="$LOG_DIR/process.log"

# Redirect stdout and stderr to the log file while still displaying it on the terminal,
# with a timestamp at the beginning of each line.
exec > >(while IFS= read -r line; do echo "$(date +'%Y-%m-%d %H:%M:%S') $line"; done | tee -a "$LOG_FILE") 2>&1

echo "Logging initiated. All output will be logged to: $LOG_FILE"

# -----------------------------------------------------------------------------
# Raw Data Setup: Make symbolic links for fastq files if not already present.
# -----------------------------------------------------------------------------

echo "Setting symbolic links for raw data in directory: $DATA"

# Check if .fastq.gz files already exist in the raw data directory.
if ls "$DATA"/*.fastq.gz 1> /dev/null 2>&1; then
    echo "Fastq files already exist in $DATA."
else
    echo "Setting symbolic links for fastq files..."
    ln -s /proj/courses/2024_fallstudie/20241209-raw_data/*.fastq.gz "$DATA"
    echo "Symbolic links set."
fi

# Report the count of files in the raw data directory.
echo "Number of files in the raw data directory: $(ls -1 "$DATA" | wc -l)"

# -----------------------------------------------------------------------------
# Script Permissions and Environment Activation
# -----------------------------------------------------------------------------

# Ensure that all shell scripts in the SCRIPTS directory are executable.
chmod +x "$SCRIPTS"/*.sh

# Activate the QIIME2 environment required for downstream analyses.
conda activate qiime2-amplicon-2024.10

# -----------------------------------------------------------------------------
# (0) Generate V3-V4 Classifier
# -----------------------------------------------------------------------------
cd "$CLASS_DIR"

#bash "$SCRIPTS/generate_classifier.sh" \


# -----------------------------------------------------------------------------
# (1) Import Amplicon Sequences into QIIME2
# -----------------------------------------------------------------------------

# Change to the import and merging subdirectory.
cd "$DATAFLOW/01-import_and_qc"

# Run the import script to bring amplicon sequences into QIIME2.
bash "$SCRIPTS/import.sh"

# -----------------------------------------------------------------------------
# (2) ASV Analysis
# -----------------------------------------------------------------------------

# Change to the ASV denoising directory.
cd "$DATAFLOW/02-ASV_denosing"

# Run the denoising script.
bash "$SCRIPTS/ASU_analysis.sh"

# Set variables for ASV quantification analysis.
ASV_RESULT_DIR_NAME="03-ASV_results" 
ASV_REP_SEQS="$DATAFLOW/02-ASV_denosing/filtered-rep-seqs.qza" 
ASV_TABLE="$DATAFLOW/02-ASV_denosing/filtered-table.qza" 

# Define parameters for ASV quantification.
ASV_SAMPLING_DEPTH=24587 
ASV_METADATA_COLUMN="sample_type"
ASV_METADATA_COLUMN2="horse"
ASV_ALPHA_RARE_MAX_DEPTH=32500
ASV_SAMPLE_TYPE="gum"
ASV_METADATA_COLUMN3="disease_state"
ASV_REF_LEVEL="healthy"
ASV_ANALYSIS_LEVEL="onset"
ASV_TAXA_COLLAPSE_LEVEL=6
ASV_TAG="ASV"

# Execute the quantification script with the defined parameters.
echo "Start ASV quantification"
bash "$SCRIPTS/run_quantification.sh" \
  "$ASV_RESULT_DIR_NAME" \
  "$ASV_REP_SEQS" \
  "$ASV_TABLE" \
  "$ASV_SAMPLING_DEPTH" \
  "$ASV_METADATA_COLUMN" \
  "$ASV_METADATA_COLUMN2" \
  "$ASV_ALPHA_RARE_MAX_DEPTH" \
  "$ASV_SAMPLE_TYPE" \
  "$ASV_METADATA_COLUMN3" \
  "$ASV_REF_LEVEL" \
  "$ASV_ANALYSIS_LEVEL" \
  "$ASV_TAXA_COLLAPSE_LEVEL" \
  "$ASV_TAG"

echo "ASV quantification completed"

# Change directory to ASV results for downstream R analysis.
cd "$DATAFLOW/03-ASV_results"   
ASV_TABLE="$DATAFLOW/02-ASV_denosing/filtered-table.qza" # Define the feature table file.

# Run the R script to perform further analysis, passing required file paths as arguments.
echo "Start ASV visualization with R"
Rscript "$SCRIPTS/run_visualization.R" "$RUNPATH/06-R" "$METADATA" "$ASV_TABLE" "$ASV_TAG"

# Remove unnecessary Rplots.pdf
if [ -f "Rplots.pdf" ]; then
    rm "Rplots.pdf"
fi

echo "ASV visualization with R completed"

# -----------------------------------------------------------------------------
# (3) OTUs Analysis
# -----------------------------------------------------------------------------

# Change to the OTUs clustering directory.
cd "$DATAFLOW/04-OTUs_clustering"

# Run the OTUs analysis script.
bash "$SCRIPTS/OTUs_analysis.sh"

# Set variables for OTUs quantification analysis.
OTUs_RESULT_DIR_NAME="05-OTUs_results"
OTUs_REP_SEQS="$DATAFLOW/04-OTUs_clustering/rep-seqs-dn-97.qza"
OTUs_TABLE="$DATAFLOW/04-OTUs_clustering/table-dn-97.qza"

# Define parameters for OTUs quantification.
OTUs_SAMPLING_DEPTH=24208
OTUs_METADATA_COLUMN="sample_type"
OTUs_METADATA_COLUMN2="horse"
OTUs_ALPHA_RARE_MAX_DEPTH=30000
OTUs_SAMPLE_TYPE="gum"
OTUs_METADATA_COLUMN3="disease_state"
OTUs_REF_LEVEL="healthy"
OTUs_ANALYSIS_LEVEL="onset"
OTUs_TAXA_COLLAPSE_LEVEL=6
OTUs_TAG="OTUs"

# Execute the OTUs quantification script with the defined parameters.
echo "Start OTUs quantification"
bash "$SCRIPTS/run_quantification.sh" \
  "$OTUs_RESULT_DIR_NAME" \
  "$OTUs_REP_SEQS" \
  "$OTUs_TABLE" \
  "$OTUs_SAMPLING_DEPTH" \
  "$OTUs_METADATA_COLUMN" \
  "$OTUs_METADATA_COLUMN2" \
  "$OTUs_ALPHA_RARE_MAX_DEPTH" \
  "$OTUs_SAMPLE_TYPE" \
  "$OTUs_METADATA_COLUMN3" \
  "$OTUs_REF_LEVEL" \
  "$OTUs_ANALYSIS_LEVEL" \
  "$OTUs_TAXA_COLLAPSE_LEVEL" \
  "$OTUs_TAG" 

echo "OTUs quantification completed"

# Change directory to ASV results for downstream R analysis.
cd "$DATAFLOW/05-OTUs_results"   

# Run the R script to perform further analysis, passing required file paths as arguments.
echo "Start OTUs visualization with R"
Rscript "$SCRIPTS/run_visualization.R" "$RUNPATH/06-R" "$METADATA" "$OTUs_TABLE" "$OTUs_TAG"

# Remove unnecessary Rplots.pdf
if [ -f "Rplots.pdf" ]; then
    rm "Rplots.pdf"
fi

echo "OTUs visualization with R completed"


# -----------------------------------------------------------------------------
# Finalization: Deactivate environment and finish processing.
# -----------------------------------------------------------------------------

conda deactivate
