#!/usr/bin/env Rscript

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
R_install <- args[1]  # R library path

# Set CRAN mirror to avoid installation errors
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Ensure R looks for packages in the specified library
.libPaths(R_install)

# Function to install missing packages
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = R_install)
  }
}

# Install required packages
install_if_missing("devtools")
install_if_missing("gplots")
install_if_missing("tidyverse")
install_if_missing("BiocManager")
install_if_missing("ggrepel")
install_if_missing("ape")

# Install Bioconductor package
BiocManager::install("ggtree")

# Load installed packages
library(devtools)
library(gplots)
library(tidyverse)
library(dplyr)
library(qiime2R)
library(ggrepel)
library(ggtree)
library(ape)

# Install qiime2R and ggrepel from GitHub
devtools::install_github("jbisanz/qiime2R")
devtools::install_github("slowkow/ggrepel")

# Load qiime2R after installation
library(qiime2R)

# Define input files
META <- args[2]  # Metadata file
TABLE <- args[3]  # Feature table (SVs)
TAG <- args[4]
TAG <- trimws(TAG)
rooted_TREE <- file.path("01-phylogeny", paste(TAG, "rooted-tree.qza", sep = "_"))
unrooted_TREE <- file.path("01-phylogeny", paste(TAG, "unrooted-tree.qza", sep = "_"))

# Define output directory
plot_dir <- file.path(getwd(), "05-visualization")
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

# Alpha Diversity Over Time ---------------------------------------------------

# Read in metadata
metadata <- read_q2metadata(META)

# Read Shannon diversity data
shannon <- read_qza("02-alpha_beta_diversity/diversity-core-metrics-phylogenetic/shannon_vector.qza")
shannon <- shannon$data %>% rownames_to_column("SampleID")

# Save Venn diagram
plot_path <- file.path(plot_dir, paste(TAG, "venn_diagram.png", sep = "_"))
png(plot_path, width = 800, height = 600)
gplots::venn(list(metadata = metadata$SampleID, shannon = shannon$SampleID))
dev.off()

# Merge metadata and diversity data
metadata <- metadata %>% left_join(shannon)

# Remove control samples
metadata <- metadata %>% filter(!(SampleID %in% c("NK", "PK")))

# Plot Shannon Diversity by Age
metadata %>%
  filter(!is.na(shannon_entropy)) %>%
  ggplot(aes(x = `age`, y = shannon_entropy, color = `sample_type`)) +
  stat_summary(geom = "errorbar", fun.data = mean_se, width = 0) +
  stat_summary(geom = "point", fun.data = mean_se) +
  xlab("Age") +
  ylab("Shannon Diversity") +
  theme_q2r() +
  scale_color_viridis_d(name = "Sample Type")
ggsave(file.path(plot_dir, paste(TAG, "Shannon_by_age.pdf", sep = "_")),
       height = 3, width = 4, device = "pdf")

# PCoA Plot -------------------------------------------------------------------
metadata <- read_q2metadata(META)
uwunifrac <- read_qza("02-alpha_beta_diversity/diversity-core-metrics-phylogenetic/unweighted_unifrac_pcoa_results.qza")
shannon <- read_qza("02-alpha_beta_diversity/diversity-core-metrics-phylogenetic/shannon_vector.qza")$data %>% rownames_to_column("SampleID")

uwunifrac$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(metadata) %>%
  left_join(shannon) %>%
  ggplot(aes(x = PC1, y = PC2, color = `sample_type`, shape = `disease_state`, size = shannon_entropy)) +
  geom_point(alpha = 0.5) +
  theme_q2r() +
  scale_shape_manual(values = c(16, 1, 17), name = "Disease state") +
  scale_size_continuous(name = "Shannon Diversity") +
  scale_color_discrete(name = "Sample type")
ggsave(file.path(plot_dir, paste(TAG, "PCoA.pdf", sep = "_")),
       height = 4, width = 5, device = "pdf")

# Taxonomy Analysis -----------------------------------------------------------
SVs <- read_qza(TABLE)$data
taxonomy_path <- file.path("03-taxonomy", paste(TAG, "taxonomy.qza", sep = "_"))
taxonomy <- read_qza(taxonomy_path)$data %>% parse_taxonomy()

taxasums <- summarize_taxa(SVs, taxonomy)$Genus

# Generate heatmap
taxa_heatmap(taxasums, metadata, "sample_type")
ggsave(file.path(plot_dir, paste(TAG, "taxa_heatmap.pdf", sep = "_")),
       height = 4, width = 8, device = "pdf")

# Generate barplot
taxa_barplot(taxasums, metadata, "sample_type")
ggsave(file.path(plot_dir, paste(TAG, "taxa_barplot.pdf", sep = "_")),
       height = 4, width = 8, device = "pdf")


# NOTE: The volcano plot and phylogenetic tree plot above can NOT be generated
#  because the qiime plug-in aldex2 is not installed in the conda environment.


# # Read differential abundance results
# results <- read_qza("04-differential_abundance/gut-test/differentials.qza")$data
# taxonomy <- read_qza("03-taxonomy/taxonomy.qza")$data
# tree <- read_qza("01-phylogeny/rooted-tree.qza")$data

# # Volcano plot
# results %>%
#   left_join(taxonomy) %>%
#   mutate(Significant = if_else(we.eBH < 0.1, TRUE, FALSE)) %>%
#   mutate(TaxonToPrint = if_else(we.eBH < 0.1, Taxon, "")) %>%
#   ggplot(aes(x = diff.btw, y = -log10(we.ep), color = Significant, label = TaxonToPrint)) +
#   geom_text_repel(size = 1, nudge_y = 0.05) +
#   geom_point(alpha = 0.6, shape = 16) +
#   theme_q2r() +
#   xlab("log2(fold change)") +
#   ylab("-log10(P-value)") +
#   theme(legend.position = "none") +
#   scale_color_manual(values = c("black", "red"))
# ggsave("05-visualization/volcano.pdf", height = 3, width = 3, device = "pdf")

# # Phylogenetic Tree Plot
# results <- results %>% mutate(Significant = if_else(we.eBH < 0.1, "*", ""))
# tree <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% results$Feature.ID])
# ggtree(tree, layout = "circular") %<+% results +
#   geom_tippoint(aes(fill = diff.btw), shape = 21, color = "grey50") +
#   geom_tiplab2(aes(label = Significant), size = 10) +
#   scale_fill_gradient2(low = "darkblue", high = "darkred", midpoint = 0, mid = "white", name = "log2(fold-change") +
#   theme(legend.position = "right")
# ggsave("05-visualization/tree.pdf", height = 10, width = 10, device = "pdf", useDingbats = FALSE)




# Import the QIIME 2 tree artifact (adjust the path as needed)
tree_artifact <- read_qza(rooted_TREE)
phylo_tree <- tree_artifact$data  # Extract the phylo object

# Plot the tree using ggtree, which utilizes ggplot2 syntax
p <- ggtree(phylo_tree) +
  geom_tiplab(size = 3, color = "darkblue") +
  ggtitle("Phylogenetic Tree from QIIME 2 Artifact") +
  theme_minimal()

# Save the plot using ggsave with a custom file name
ggsave(file.path(plot_dir, paste(TAG, "rooted-tree.pdf", sep = "_")),
       plot = p,
       height = 4, width = 8, device = "pdf")


tree_artifact <- read_qza(unrooted_TREE)
phylo_tree <- tree_artifact$data  # Extract the phylo object

# Plot the tree using ggtree, which utilizes ggplot2 syntax
p <- ggtree(phylo_tree) +
  geom_tiplab(size = 3, color = "darkblue") +
  ggtitle("Phylogenetic Tree from QIIME 2 Artifact") +
  theme_minimal()

# Save the plot using ggsave with a custom file name
ggsave(file.path(plot_dir, paste(TAG, "unrooted-tree.pdf", sep = "_")),
       plot = p,
       height = 4, width = 8, device = "pdf")

