# Amplicon Analysis Pipeline with QIIME 2: ASV and OTUs Approaches

## Overview

This pipeline performs amplicon sequence analysis using **QIIME 2**, supporting both **OTUs (Operational Taxonomic Unit)** clustering 
and **ASV (Amplicon Sequence Variant)** inference. It includes all steps from primer trimming to taxonomic assignment and integrates downstream visualization using 
- **R** (e.g.,qiime2R, ggplot2`).
- **Python**
- and the online visualization of qiime "https://view.qiime2.org/"
---

## Features

- ASV inference using `dada2`
- Trains own classifier for the V3-V4 region of the 16S rRNA gene
- Sample and feature filtering based on metadata
- OTUs clustering at 97% identity using `vsearch` de-novo method
- Visualization and statistical analysis in R using `qiime2R`

---

## Requirements

### Software

- [QIIME 2 (2024.x or later)](https://qiime2.org/)
- [R (>= 4.0)](https://cran.r-project.org/)
- R packages: `qiime2R`, `ggplot2`, `dplyr`, `devtools`, `ggrepel`, etc.
- Python: `os`, `numpy`, `pandas`, etc. 
- Bash


### Data

- Demultiplexed paired-end FASTQ files
- Metadata file in TSV format (`metadata.tsv`)

---

## Running the script

To run the script, first navigate to the directory where run.sh is located and run the command:

bash -i run.sh

Ensure that QIIME 2 is installed in a Conda environment named qiime2-amplicon-2024.10, including all necessary plugins such as DADA2 and ALDEx2.
If you havent installed it, a yml-file can be found in "05-conda_env" to export the qiime2 Conda environment. Also the needed R packages are
stored under "06-R", if you havent install them in your regular directory. 

The classifier for the taxonomic analysis of the V3-V4 region will be generted at first by the script which will take a couple hours, so be patient. 

You can either place your FASTQ files in the local directory 01-data/20241009-raw_data, or—if you're working on the FHWN server and have access to 
/proj/courses/2024_fallstudie/20241209-raw_data/—you can simply run the script, and it will automatically mirror the files to the appropriate local directory.

The pipeline is available on GitHub and can be cloned with the command:

`git clone https://github.com/besAvdiu303/EORTH_qiime2_analysis.git`

