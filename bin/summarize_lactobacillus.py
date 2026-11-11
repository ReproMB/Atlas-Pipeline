#!/usr/bin/env python3
"""
summarize_lactobacillus.py

Reads QIIME2 feature table TSV, taxonomy TSV, and metadata TSV.
Filters to Lactobacillus ASVs, calculates relative abundance per sample.
Generates a CSV summary and simple barplot/boxplot.

Usage:
    summarize_lactobacillus.py --table feature_table.tsv --tax taxonomy.tsv \
        --meta metadata.tsv --prefix study_id
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--table', required=True)
    parser.add_argument('--tax', required=True)
    parser.add_argument('--meta', required=True)
    parser.add_argument('--prefix', required=True)
    args = parser.parse_args()

    # Load table (skip QIIME2 header rows)
    # The first line is usually # Constructed from biom file, second is #OTU ID...
    with open(args.table) as f:
        skip_rows = 1 if f.readline().startswith('# Constructed') else 0
        
    df_table = pd.read_csv(args.table, sep='\t', skiprows=skip_rows)
    df_table.rename(columns={df_table.columns[0]: 'Feature ID'}, inplace=True)
    df_table.set_index('Feature ID', inplace=True)
    
    # Load Taxonomy
    df_tax = pd.read_csv(args.tax, sep='\t')
    
    # Load Metadata
    df_meta = pd.read_csv(args.meta, sep='\t')
    if 'sample-id' not in df_meta.columns:
        df_meta.rename(columns={df_meta.columns[0]: 'sample-id'}, inplace=True)
        
    # Check if the table is already collapsed by taxonomy
    # If the intersection with taxonomy feature IDs is empty, the index is the taxonomy string
    intersection = df_table.index.intersection(df_tax['Feature ID'])
    
    if len(intersection) > 0:
        # Uncollapsed table (ASVs) — filter to d__Bacteria features only
        bacteria_features = df_tax[
            df_tax['Taxon'].str.startswith('d__Bacteria', na=False)
        ]['Feature ID'].tolist()
        df_table = df_table.loc[df_table.index.intersection(bacteria_features)]

        lacto_features = df_tax[
            df_tax['Taxon'].str.contains('Lactobacillus', case=False, na=False) &
            df_tax['Taxon'].str.startswith('d__Bacteria', na=False)
        ]['Feature ID'].tolist()
        df_lacto = df_table.loc[df_table.index.intersection(lacto_features)]
    else:
        # Collapsed table (Taxonomy strings) — keep only d__Bacteria rows
        df_table = df_table[df_table.index.str.startswith('d__Bacteria', na=False)]
        df_lacto = df_table[df_table.index.str.contains('Lactobacillus', case=False, na=False)]
        
    # Calculate relative abundance (denominator is bacteria-only counts)
    # 1. Total bacteria counts per sample
    total_counts = df_table.sum(axis=0)
    
    # 2. Lacto counts per sample
    lacto_counts = df_lacto.sum(axis=0)
    
    # 3. Relative abundance
    rel_abund = (lacto_counts / total_counts).fillna(0).reset_index()
    rel_abund.columns = ['sample-id', 'Lactobacillus_RelAbund']
    
    # Merge with metadata
    df_final = pd.merge(rel_abund, df_meta, on='sample-id', how='inner')
    
    # Save CSV
    out_csv = f"{args.prefix}_lactobacillus_summary.csv"
    df_final.to_csv(out_csv, index=False)
    print(f"Saved {out_csv}")
    
    # Generate Boxplot
    if 'ConditionAtlas' in df_final.columns and len(df_final) > 0:
        plt.figure(figsize=(10, 6))
        sns.boxplot(data=df_final, x='ConditionAtlas', y='Lactobacillus_RelAbund', hue='SiteAtlas')
        sns.stripplot(data=df_final, x='ConditionAtlas', y='Lactobacillus_RelAbund', hue='SiteAtlas', dodge=True, color='black', alpha=0.5)
        plt.title(f"{args.prefix} - Lactobacillus Relative Abundance")
        plt.ylabel("Relative Abundance")
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig(f"{args.prefix}_lactobacillus_boxplot.png", dpi=300)
        print(f"Saved {args.prefix}_lactobacillus_boxplot.png")

if __name__ == "__main__":
    main()
