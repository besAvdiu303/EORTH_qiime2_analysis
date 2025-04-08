#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ----------------------------
# Parse command-line arguments
# ----------------------------
parser = argparse.ArgumentParser(description="Generate and save genus clustermap.")
parser.add_argument("tsv_path", type=str, help="Path to the input TSV file")
parser.add_argument("output_path", type=str, help="Path to save the PDF (including filename)")
args = parser.parse_args()

# Load and preprocess data
abundances = pd.read_table(args.tsv_path, skiprows=1, index_col=0)
abundances.index = abundances.index.str.split(";").str[5]
abundances = abundances[~abundances.index.isin(["g__", "__"])]
abundances = abundances.sample(50)


transformed = abundances.apply(
    lambda xs: np.log(xs + 0.5) -np.log(xs.mean() + 0.5),
    axis=0)


# ----------------------------
# Plot and save
# ----------------------------
clustergrid = sns.clustermap(transformed.T, cmap="magma", xticklabels=True, figsize=(18, 6))

# Ensure output directory exists
os.makedirs(os.path.dirname(args.output_path), exist_ok=True)

# Save the figure
clustergrid.savefig(args.output_path)
plt.close()

print(f"Clustermap saved to: {args.output_path}")