# Amplicon Analysis Pipeline with QIIME 2: ASV and OTUs Approaches

## Overview

This pipeline performs amplicon sequence analysis using **QIIME 2**, supporting both **OTUs (Operational Taxonomic Unit)** clustering 
and **ASV (Amplicon Sequence Variant)** inference. It includes all steps from primer trimming to taxonomic assignment and integrates downstream visualization using **R** (e.g., ggplot2`).

---

## Features

- ASV inference using DADA2
- Sample and feature filtering based on metadata
- Primer removal using `cutadapt`
- Paired-end merging using `vsearch`
- Quality filtering via `q-score`
- Dereplication and contaminant removal
- OTUs clustering at 97% identity
- Sample and feature filtering based on metadata
- Visualization and statistical analysis in R

---

## Requirements

### Software

- [QIIME 2 (2024.x or later)](https://qiime2.org/)
- [R (>= 4.0)](https://cran.r-project.org/)
- R packages: `qiime2R`, `ggplot2`, `dplyr`, `devtools`, `ggrepel`, etc.
- Bash


### Data

- Demultiplexed paired-end FASTQ files
- Metadata file in TSV format (`metadata.tsv`)

---

## Running the script

To run the script, first navigate to the directory where run.sh is located.
Ensure that QIIME 2 is installed in a Conda environment named qiime2-amplicon-2024.10, including all necessary plugins such as DADA2 and ALDEx2.
If you havent installed it, a yml-file can be found in 05-conda_env to export the qiime2 Conda environment. Also the needed R packages are
stored under 06-R, if you havent install them in your regular directory. 

You can either place your FASTQ files in the local directory 01-data/20241009-raw_data, or—if you're working on the FH server and have access to 
/proj/courses/2024_fallstudie/20241209-raw_data/—you can simply run the script, and it will automatically mirror the files to the appropriate local directory.

The pipeline can be clone from GitHub with the command:

git clone https://github.com/besAvdiu303/EORTH_qiime2_analysis.git

