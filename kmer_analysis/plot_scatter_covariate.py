"""
Author: Julie Chen
plot_scatter_covariate.py

Given a user-specified csv file containing enrichment scores for motifs of interest from extract_motifs.py
and a covariate column of interest, generates scatter plots plotting covariate against enrichment of motifs
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import os
from adjustText import adjust_text

import numpy as np
from scipy.stats import spearmanr
import seaborn as sns

def shorten_species_name(name):
    """Shorten species name for plotting, e.g. Archaeoglobus_fulgidus -> A. fulgidus"""
    parts = name.split('_')
    if len(parts) >= 2:
        return f"{parts[0][0]}. {parts[1]}"
    else:
        return name

def get_motif_columns(df):
    """
    Return a list of motif columns starting from 'tRNA_scan_motif',
    excluding any columns ending in '_sequence' or '_kmer'.
    """
    start_idx = df.columns.get_loc('tRNA_scan_motif')
    motif_cols = []
    for col in df.columns[start_idx:]:
        if col.endswith('_sequence') or col.endswith('_kmer'):
            continue
        motif_cols.append(col)
    return motif_cols

def plot_covariate_vs_motif(df, covariate, motif_col, output_prefix):
    """
    Make a scatter plot of covariate vs motif enrichment,
    labeling points with shortened species names.
    """
    #plt.rcParams['font.family'] = 'Helvetica'
    plt.figure(figsize=(7.5,6))
    plt.scatter(df[covariate], df[motif_col], color="#5473CF")

    texts = []
    for _, row in df.iterrows():
        if pd.notna(row[motif_col]):
            label = shorten_species_name(row['species'])
            texts.append(plt.text(row[covariate], row[motif_col], label,
                                fontsize=9, ha='center', va='bottom', zorder=2, fontstyle='italic')) 

    # Adjust labels to avoid overlaps, allow both x and y movement
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='lightgray', alpha=0.5, lw=0.7), 
                expand=(1.25, 1.25))
    
    ax = plt.gca()  # get current axes

    # Scatter points + linear regression with 95% CI (funnel)
    sns.regplot(
        x=covariate, y=motif_col, data=df,
        scatter=True, ci=95, color="#5473CF",
        line_kws={'color':"#0B2779", 'linestyle':':', 'linewidth':2},
        scatter_kws={'zorder':2}
    )

    # Hide top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Spearman correlation
    x = df[covariate].values
    y = df[motif_col].values
    mask = ~np.isnan(x) & ~np.isnan(y)
    rho, pval = spearmanr(x[mask], y[mask])
    # Annotate Spearman rho and p-value
    plt.text(0.95, 0.05, f"$\\rho$ = {rho:.2f}\np = {pval:.3g}",
            transform=ax.transAxes, fontsize=12, fontname="Helvetica", verticalalignment='bottom', horizontalalignment='right')
    
    if covariate == "nap" or "nap_noFerritin":
        cov_name = f"NAP abundance (%)"
    elif covariate == "nap_known" or "nap_known_noFerritin":
        cov_name = f"Known NAP abundance (%)"
    else:
        cov_name = covariate

    plt.xlabel(cov_name, fontsize=14, fontname="Helvetica")
    plt.ylabel(r'$\log_{2}$' + '(real / randomized)', fontsize=14, fontname="Helvetica")

    plt.title(f"{motif_col}", fontsize=22, fontname="Helvetica")
    plt.axhline(0, color='gray', linestyle='--', linewidth=2)
    plt.tick_params(axis='x', labelsize=12, labelfontfamily="Helvetica")
    plt.tick_params(axis='y', labelsize=12, labelfontfamily="Helvetica")
    plt.tight_layout()

    output_file = f"{output_prefix}_{motif_col}_vs_{covariate}.png"
    plt.savefig(output_file, dpi=1200)
    plt.close()
    print(f"Saved plot: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Scatter plots of covariate vs motif enrichments")
    parser.add_argument('csv_file', type=str, help="CSV file containing enrichment summary")
    parser.add_argument('covariate', type=str, help="Covariate column name to plot against")
    args = parser.parse_args()

    # Load CSV
    df = pd.read_csv(args.csv_file)
    if 'species' not in df.columns:
        raise ValueError("CSV must contain a 'species' column")

    # Identify motif columns dynamically
    motif_cols = get_motif_columns(df)

    # Prepare output prefix from csv filename
    base_name = os.path.basename(args.csv_file)
    output_prefix = os.path.splitext(base_name)[0]

    # Generate scatter plots for each motif
    for motif_col in motif_cols:
        if motif_col not in df.columns:
            print(f"Warning: motif column {motif_col} not found, skipping")
            continue
        plot_covariate_vs_motif(df, args.covariate, motif_col, output_prefix)

if __name__ == "__main__":
    main()
