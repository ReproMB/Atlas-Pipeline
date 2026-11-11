#!/usr/bin/env python3
"""
parse_study_metadata.py

Finds and parses the metadata file for a specific study inside ESTUDIOS/.
Filters samples that have both Condition and Site.
Extracts Age, BMI, Ethnicity and Experiment (if present).
Generates metadata.tsv for QIIME2 and sample_filter.txt for the manifest.

Usage:
    parse_study_metadata.py --study PRJNA329174 --metadata-dir ../ESTUDIOS
"""

import argparse
import pandas as pd
import numpy as np
import os
import glob
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--study', required=True)
    parser.add_argument('--metadata-dir', required=True)
    parser.add_argument('--out-tsv', default='metadata.tsv')
    parser.add_argument('--out-filter', default='sample_filter.txt')
    parser.add_argument('--sites', default=None,
                        help='Comma-separated body sites for this processing unit. '
                             'When given, only samples whose SiteAtlas matches one of '
                             'them are kept, so a multi-site study does not put the '
                             'same samples into every site unit.')
    args = parser.parse_args()

    # Find the metadata file (it might have spaces or be .csv/.xlsx)
    study_id = args.study.strip()
    search_pattern = os.path.join(args.metadata_dir, f"*{study_id}*.*")
    candidates = glob.glob(search_pattern)
    
    valid_files = [f for f in candidates if f.endswith('.xlsx') or f.endswith('.csv') or f.endswith('.tsv')]
    
    if not valid_files:
        print(f"WARNING: No metadata file found for {study_id} in {args.metadata_dir}", file=sys.stderr)
        # Create empty files to avoid pipeline crashing
        with open(args.out_tsv, 'w') as f:
            f.write("sample-id\n")
        with open(args.out_filter, 'w') as f:
            pass
        sys.exit(0)
        
    metadata_file = valid_files[0]
    print(f"Found metadata file: {metadata_file}", file=sys.stderr)
    
    if metadata_file.endswith('.csv'):
        # Auto-detect separator: some studies use ; instead of ,
        with open(metadata_file) as f:
            first_line = f.readline()
        sep = ';' if first_line.count(';') > first_line.count(',') else ','
        df = pd.read_csv(metadata_file, sep=sep)
    elif metadata_file.endswith('.tsv'):
        df = pd.read_csv(metadata_file, sep='\t')
    else:
        df = pd.read_excel(metadata_file)
        
    # Standardize column names (lowercase for matching)
    col_map = {c.lower().strip(): c for c in df.columns}
    
    # Identify key columns
    id_col = col_map.get('run') or col_map.get('sample_id') or col_map.get('sample name') or col_map.get('sample id') or col_map.get('experiment')
    condition_col = col_map.get('conditionatlas') or col_map.get('condition')
    site_col = col_map.get('siteatlas') or col_map.get('site')
    age_col = col_map.get('ageatlas') or col_map.get('age')
    bmi_col = col_map.get('bmiatlas') or col_map.get('bmi')
    ethnicity_col = col_map.get('ethnicityatlas') or col_map.get('ethnicity')
    # Menstrual cycle phase — matched BY NAME (like every other column), never
    # by position. Tolerates the legacy 'MentrualCycle' typo just in case.
    menstrual_col = col_map.get('menstrualcycle') or col_map.get('mentrualcycle')
    experiment_col = col_map.get('experiment')
    
    if not id_col or not condition_col or not site_col:
        print(f"WARNING: Missing required columns (ID, Condition, Site) in {metadata_file}", file=sys.stderr)
        with open(args.out_tsv, 'w') as f:
            f.write("sample-id\n")
        with open(args.out_filter, 'w') as f:
            pass
        sys.exit(0)

    # Filter out rows missing condition or site
    df_filtered = df.dropna(subset=[condition_col, site_col])

    # Keep only the sites belonging to this processing unit. Without this a
    # multi-site study hands every unit the full sample set, so the same
    # samples end up counted once per site.
    if args.sites:
        wanted = {s.strip().lower() for s in args.sites.split(',') if s.strip()}
        if wanted:
            site_norm = df_filtered[site_col].astype(str).str.strip().str.lower()
            before = len(df_filtered)
            df_filtered = df_filtered[site_norm.isin(wanted)]
            print(f"Site filter {sorted(wanted)}: {before} -> {len(df_filtered)} samples",
                  file=sys.stderr)

    out_data = []
    filter_ids = []
    
    for _, row in df_filtered.iterrows():
        sample_id = str(row[id_col]).strip()
        condition = str(row[condition_col]).strip()
        site = str(row[site_col]).strip()
        
        # Parse optional columns
        age = str(row[age_col]).strip() if age_col and pd.notna(row[age_col]) else "NA"
        bmi = str(row[bmi_col]).strip() if bmi_col and pd.notna(row[bmi_col]) else "NA"
        ethnicity = str(row[ethnicity_col]).strip() if ethnicity_col and pd.notna(row[ethnicity_col]) else "NA"
        menstrual = str(row[menstrual_col]).strip() if menstrual_col and pd.notna(row[menstrual_col]) else "NA"
        experiment = str(row[experiment_col]).strip() if experiment_col and pd.notna(row[experiment_col]) else "NA"

        # Save to list
        out_data.append({
            'sample-id': sample_id,
            'ConditionAtlas': condition,
            'SiteAtlas': site,
            'AgeAtlas': age,
            'BMIAtlas': bmi,
            'EthnicityAtlas': ethnicity,
            'MenstrualCycle': menstrual,
            'Experiment': experiment
        })
        
        filter_ids.append(sample_id)
        if experiment != "NA" and experiment != sample_id:
            filter_ids.append(experiment) # Allow filtering by experiment ID as well
        # For non-SRA studies (e.g. Chinese databases like OMIX/GSA),
        # filenames often use sample_name or run_name instead of sample_id.
        # Add these as alternative filter IDs so the manifest can match them.
        for alt_col_name in ['run_name', 'sample_name', 'run_id']:
            alt_col_real = col_map.get(alt_col_name.lower())
            if alt_col_real and pd.notna(row.get(alt_col_real)):
                alt_val = str(row[alt_col_real]).strip()
                if alt_val and alt_val != sample_id and alt_val != "NA":
                    filter_ids.append(alt_val)
            
    out_df = pd.DataFrame(out_data)
    
    if out_df.empty:
        print(f"WARNING: No samples passed the filter in {metadata_file}", file=sys.stderr)
        with open(args.out_tsv, 'w') as f:
            f.write("sample-id\n")
        with open(args.out_filter, 'w') as f:
            pass
        sys.exit(0)
        
    # Write QIIME2 metadata TSV
    out_df.to_csv(args.out_tsv, sep='\t', index=False)
    
    # Write sample filter
    with open(args.out_filter, 'w') as f:
        for sid in set(filter_ids):
            f.write(f"{sid}\n")
            
    print(f"Extracted {len(out_df)} valid samples to {args.out_tsv}", file=sys.stderr)

if __name__ == "__main__":
    main()
