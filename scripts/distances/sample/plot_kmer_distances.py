#!/usr/bin/env python3
"""
plot_kmer_distances.py

Enhanced plotting for distance matrices and PCoA coordinates, with:
- auto-stretched color scale for each heatmap to show subtle differences
- consistent sample ordering using hierarchical clustering
- PCoA axis labels include % explained variance if provided
- region and grouping names added to PCoA titles and filenames
"""

import argparse
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

def plot_distance_matrix(dm, title, out_path, vmin=None, vmax=None, figsize=(10,8), cmap="viridis"):
    plt.figure(figsize=figsize)
    ax = sns.heatmap(
        dm, square=True, cmap=cmap,
        cbar_kws={'label':'Distance'},
        vmin=vmin, vmax=vmax,
        annot=False
    )
    
    # Automatically scale font size based on number of samples
    n_samples = len(dm.columns)
    font_size = max(4, min(10, 200 // n_samples))  # keeps it between 4 and 10
    
    # Set all x and y tick labels
    ax.set_xticks(np.arange(n_samples) + 0.5)
    ax.set_xticklabels(dm.columns, rotation=90, fontsize=font_size)
    ax.set_yticks(np.arange(n_samples) + 0.5)
    ax.set_yticklabels(dm.index, rotation=0, fontsize=font_size)
    
    plt.title(title)
    plt.tight_layout()
    # Reverse the y-axis so smallest y_position appears at the bottom
    ax.invert_yaxis()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Distance heatmap saved to {out_path}")

def plot_pcoa(coords, eigenvals=None, grouping=None, title="PCoA plot", out_path="pcoa_plot.pdf", figsize=(8,6)):
    plt.figure(figsize=figsize)
    sns.set(style="whitegrid")

    if grouping is not None:
        coords = coords.merge(grouping, left_index=True, right_on='sample')
        hue_col = 'group'
    else:
        hue_col = None

    if hue_col:
        ax = sns.scatterplot(
            data=coords,
            x=coords.columns[0],
            y=coords.columns[1],
            hue=hue_col,
            s=100,
            palette="tab10"
        )
    else:
        ax = sns.scatterplot(
            data=coords,
            x=coords.columns[0],
            y=coords.columns[1],
            s=100
        )

    for i, row in coords.iterrows():
        ax.text(
            row[coords.columns[0]] + 0.01,
            row[coords.columns[1]] + 0.01,
            str(row.get('sample', i)),
            fontsize=6,
            alpha=0.8
        )

    if eigenvals is not None:
        total_var = sum(eigenvals)
        explained = [(100*ev/total_var) for ev in eigenvals[:2]]
        ax.set_xlabel(f"{coords.columns[0]} ({explained[0]:.1f}%)")
        ax.set_ylabel(f"{coords.columns[1]} ({explained[1]:.1f}%)")
    else:
        ax.set_xlabel(coords.columns[0])
        ax.set_ylabel(coords.columns[1])

    ax.set_title(title)
    if hue_col:
        plt.legend(title=hue_col, bbox_to_anchor=(1.05,1), loc='upper left')
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"PCoA plot saved to {out_path}")

def main(args):
    os.makedirs(args.output_dir, exist_ok=True)

    # region is user-specified or defaults to TprK_V_regions
    region = args.region if args.region else 'TprK_V_regions'

    # --- NEW: load y-position order ---
    ypos_file = args.ypos
    ypos_df = pd.read_csv(ypos_file, dtype=str)
    ypos_df['y_position'] = pd.to_numeric(ypos_df['y_position'], errors='coerce')
    ypos_df = ypos_df.sort_values('y_position')
    sample_order = ypos_df['sample'].tolist()
    print(f"Loaded y_position order from {ypos_file}: {len(sample_order)} samples")

    # Load all distance matrices
    all_matrices = []
    for dist_file in args.dist_files:
        dm = pd.read_csv(dist_file, sep="\t", index_col=0)
        all_matrices.append((dist_file, dm))

    # --- REPLACE cluster order with fixed y_position order ---
    for dist_file, dm in all_matrices:
        # filter only samples that exist in dm
        ordered = [s for s in sample_order if s in dm.index]
        missing = set(dm.index) - set(ordered)
        if missing:
            print(f"Warning: samples missing from y_position file: {', '.join(missing)}")

        dm_ordered = dm.loc[ordered, ordered]

        # dynamic vmin/vmax scaling
        vals = dm_ordered.values[np.triu_indices_from(dm_ordered.values, k=1)]
        vmin, vmax = vals.min(), vals.max()
        margin = (vmax - vmin) * 0.01
        vmin -= margin
        vmax += margin

        base_name = os.path.basename(dist_file).replace(".tsv", "")
        title = f"{region} {base_name}"
        out_name = f"{region}_{base_name}_heatmap.pdf"
        out_path = os.path.join(args.output_dir, out_name)
        plot_distance_matrix(dm_ordered, title, out_path, vmin=vmin, vmax=vmax)

    # Optional eigenvalues
    eigenvals = None
    if args.eigenvalues_csv:
        eigenvals = pd.read_csv(args.eigenvalues_csv, header=None).iloc[:,0]
        eigenvals = pd.to_numeric(eigenvals, errors='coerce').dropna().tolist()

    # Plot PCoA for no grouping
    coords = pd.read_csv(args.pcoa_csv, index_col=0)
    pcoa_title = f"{region} PCoA"
    pcoa_out = os.path.join(args.output_dir, f"{region}_pcoa_plot.pdf")
    plot_pcoa(coords, eigenvals=eigenvals, grouping=None, title=pcoa_title, out_path=pcoa_out)

    # Plot PCoA for each grouping file
    if args.grouping_files:
        for gfile in args.grouping_files:
            grouping_base = os.path.basename(gfile).replace(".csv","")
            grouping = pd.read_csv(gfile, dtype=str)
            if 'sample' not in grouping.columns or 'group' not in grouping.columns:
                raise ValueError(f"Grouping CSV {gfile} must have columns: sample, group")
            pcoa_title = f"{region} {grouping_base} PCoA"
            pcoa_out = os.path.join(args.output_dir, f"{region}_{grouping_base}_pcoa_plot.pdf")
            plot_pcoa(coords, eigenvals=eigenvals, grouping=grouping, title=pcoa_title, out_path=pcoa_out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dist_files", nargs='+', required=True, help="TSV distance matrices (Bray, JSD, CLR, etc.)")
    parser.add_argument("--pcoa_csv", required=True, help="PCoA coordinates CSV")
    parser.add_argument("--grouping_files", nargs='*', help="Optional list of CSVs with sample,group for coloring")
    parser.add_argument("--eigenvalues_csv", help="Optional CSV with PCoA eigenvalues for explained variance")
    parser.add_argument("--output_dir", default="plots", help="Directory to save plots")
    parser.add_argument("--region", help="Optional region name to override default ('TprK_V_regions')")
    parser.add_argument("--ypos", help="Optional ypositions file to override default ('ypositions.csv')")
    args = parser.parse_args()
    main(args)
