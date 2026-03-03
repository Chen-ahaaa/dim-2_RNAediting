#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  8 14:38:38 2025

@author: cxn
""" 
import os
import glob
import pandas as pd
from collections import defaultdict
import argparse

############################################################
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_dir", required=True)
parser.add_argument("-o", "--output_dir", required=True)
parser.add_argument("--sep", default="\t")
args = parser.parse_args()

def merge_vs_groups(input_dir, output_dir, sep="\t"):
    os.makedirs(output_dir, exist_ok=True)

    all_files = glob.glob(os.path.join(input_dir, "*_vs_*.edsite_ann.csv"))
    if not all_files:
        print("No _vs_ files found.")
        return

    prefix_groups = defaultdict(list)
    for f in all_files:
        fname = os.path.basename(f)
        prefix = fname.split("_vs_")[0]
        prefix_groups[prefix].append(f)

    for prefix, files in prefix_groups.items():
        print(f"Processing prefix: {prefix} ({len(files)} files)")
        dfs = [pd.read_csv(f, sep=sep, dtype=str).fillna("") for f in files]
        if not dfs:
            continue
        merged_df = pd.concat(dfs, ignore_index=True)

        all_cols = merged_df.columns.tolist()
        agg_cols = [c for c in all_cols if (c.startswith("g") or c == "DNA_sample") and "gene" not in c]
        group_cols = [c for c in all_cols if c not in agg_cols]                    
        
        agg_dict = {c: lambda x: ";".join([v for v in x if v != ""]) for c in agg_cols}
        merged_df = (
            merged_df
            .groupby(group_cols, dropna=False, as_index=False)
            .agg(agg_dict))
        
        for c in agg_cols:
            merged_df[c] = merged_df[c].str.strip(";")
            merged_df = merged_df[merged_df["DNA_sample"].str.split(";").apply(len) == len(dfs)]
    
        out_path = os.path.join(output_dir, f"{prefix}_merged.edsite_ann.csv")
        merged_df.to_csv(out_path, sep=sep, index=False)

def merge_by_sample_key(input_dir, output_dir, sep="\t"):
    os.makedirs(output_dir, exist_ok=True)

    # find all *_merged.txt files
    all_files = glob.glob(os.path.join(input_dir, "*_merged.edsite_ann.csv"))
    if not all_files:
        print("No *_merged.edsite_ann.csv files found.")
        return

    # groupby *_merged.txt 
    prefix_groups = defaultdict(list)
    for f in all_files:
        fname = os.path.basename(f)
        base = fname.rsplit("_merged.edsite_ann.csv", 1)[0]
        if "_" in base:
            key = "_".join(base.split("_")[1:])  # '6dpf_as_1'
        else:
            key = base
        prefix_groups[key].append(f)

    for key, files in prefix_groups.items():
        print(f"Processing group: {key} ({len(files)} files)")

        merged_df = None
        for idx, f in enumerate(files):
            df = pd.read_csv(f, sep=sep, dtype=str)
            if idx == 0:
                merged_df = df
            else:
                merged_df = pd.concat([merged_df, df], ignore_index=True)

        if merged_df is not None:
            out_path = os.path.join(output_dir, f"{key}_final.edsite_ann.csv")
            merged_df.to_csv(out_path, sep=sep, index=False)
            
merge_vs_groups(args.input_dir,args.input_dir)
merge_by_sample_key(args.input_dir, args.output_dir)
