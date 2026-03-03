#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 18 09:48:01 2025

@author: cxn
"""
import pandas as pd
import os
import subprocess
import pysam
import argparse
import tempfile
import ast
import pyranges as pr
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input",type=str,required=True)
parser.add_argument("-o","--outdir",type=str,required=True)
parser.add_argument("-d","--strand_source",type=int,required=True,help="Use 0 for sense strand and 1 for antisense strand.")
parser.add_argument("-s","--spacies",type=str,required=True)
parser.add_argument("-n","--sample_name",type=str,required=True)
parser.add_argument("-g","--gff",type=str,required=True)
parser.add_argument("-D","--DNA_name",type=str,default="-")
parser.add_argument("-fq","--frequency",type=float,default=0)
parser.add_argument("-r","--region",type=str,default=None)
parser.add_argument("-header","--vcf_header",type=str,default=None)


args = parser.parse_args()
editing_file=args.input
out_dir=args.outdir
snpeff_spacies=args.spacies
strand_source=args.strand_source
sample_name=args.sample_name
DNA_name=args.DNA_name
gff_file=args.gff
min_freq=args.frequency
region=args.region
vcf_header=args.vcf_header
############################################################################################
def filter_ed_by_gff(ed_df, gff_df):
    gff_pr = pr.PyRanges(
        pd.DataFrame({
            "Chromosome": gff_df[0],
            "Start": gff_df[3].astype(int),
            "End": gff_df[4].astype(int)}))
    ed_range = ed_df.assign(
        Chromosome=ed_df["Region"],
        Start=ed_df["Position"].astype(int),
        End=ed_df["Position"].astype(int) + 1)
    ed_pr = pr.PyRanges(ed_range)
    intersected = ed_pr.join(gff_pr)
    if not intersected.df.empty:
        filtered_rows = intersected.df.loc[:, ed_df.columns]
    else:
        filtered_rows = ed_df.iloc[0:0]
    filtered_rows = filtered_rows.drop_duplicates()
    return filtered_rows
############################################################################
ed_df_all=pd.read_csv(editing_file,dtype=str,sep="\t")
gff_df=pd.read_csv(gff_file,comment="#",sep="\t",dtype=str,header=None)
gff_df=gff_df[gff_df[2]!="region"]
ed_df=ed_df_all[(ed_df_all["Frequency"].astype(float)>min_freq)]
mask = ed_df.apply(lambda r: str(r["AllSubs"]) in str(r["gAllSubs"]), axis=1)
ed_df = ed_df.loc[~mask].copy()

if region:
    region_pr = pr.read_bed(region)
    ed_range = ed_df.assign(Chromosome=ed_df["Region"],Start=ed_df["Position"].astype(int),End=ed_df["Position"].astype(int) + 1)
    ed_pr = pr.PyRanges(ed_range)
    intersected = ed_pr.join(region_pr)
    ed_df = (intersected.df.loc[:, ed_df.columns] if not intersected.df.empty else ed_df.iloc[0:0])
    
if strand_source==0:
    ed_df = ed_df[ed_df["AllSubs"] == "AG"]
    ed_df = ed_df.assign(sup_reads = ed_df["BaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[2]))
    ed_df = ed_df.assign(gsup_reads = ed_df["gBaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[2]))
    gff_df=gff_df[gff_df[6]=="+"]
    filtered_rows=filter_ed_by_gff(ed_df, gff_df)
    filtered_rows["strand"]="+"
elif strand_source==1:
    ed_df = ed_df[ed_df["AllSubs"] == "TC"]
    ed_df = ed_df.assign(sup_reads = ed_df["BaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[1]))
    ed_df = ed_df.assign(gsup_reads = ed_df["gBaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[1]))
    gff_df=gff_df[gff_df[6]=="-"]
    filtered_rows=filter_ed_by_gff(ed_df, gff_df)
    filtered_rows["strand"]="-"
elif strand_source==2:
    ed_df_sense = ed_df[ed_df["AllSubs"] == "AG"]
    ed_df_sense = ed_df_sense.assign(sup_reads = ed_df_sense["BaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[2]))
    ed_df_sense = ed_df_sense.assign(gsup_reads = ed_df_sense["gBaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[2]))
    gff_df_sense=gff_df[gff_df[6]=="+"]
    filtered_sense=filter_ed_by_gff(ed_df_sense, gff_df_sense)
    filtered_sense["strand"]="+"
    
    ed_df_antisense = ed_df[ed_df["AllSubs"] == "TC"]
    ed_df_antisense = ed_df_antisense.assign(sup_reads = ed_df_antisense["BaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[1]))
    ed_df_antisense = ed_df_antisense.assign(gsup_reads = ed_df_antisense["gBaseCount[A,C,G,T]"].apply(lambda x: "-" if x == "-" else ast.literal_eval(x)[1]))
    gff_df_antisense=gff_df[gff_df[6]=="-"]
    filtered_antisense=filter_ed_by_gff(ed_df_antisense, gff_df_antisense)
    filtered_antisense["strand"]="-"
    filtered_rows = pd.concat([filtered_sense, filtered_antisense], ignore_index=True)


annot_list = []
with tempfile.TemporaryDirectory() as tmpdir:
    vcf_path = os.path.join(tmpdir, f"{sample_name}.vcf")
    ann_vcf_path = os.path.join(tmpdir, f"{sample_name}.ann.vcf")
    genes_stats_path = os.path.join(tmpdir, "snpEff_genes.txt")
    html_stats_path = os.path.join(tmpdir, "snpEff_summary.html")
    
    # construct tmp vcf
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        if vcf_header and os.path.exists(vcf_header):
            with open(vcf_header) as header_file:
                for line in header_file:
                    line = line.rstrip()
                    f.write(line + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for idx, row in filtered_rows.iterrows():
            f.write(f"{row['Region']}\t{row['Position']}\t.\t{row['Reference']}\t{row['AllSubs'][1]}\t.\t.\t.\n")
    # get ann.vcf
    subprocess.run([
        "java", "-jar", "/mnt/san3/usr/cxn/soft/snpEff/snpEff.jar",
        "-ud", "0",
        "-csvStats", genes_stats_path,
        "-stats", html_stats_path,
        "-v", snpeff_spacies,
        vcf_path
    ], stdout=open(ann_vcf_path, "w"), check=True)
    
    # 3. 解析注释 VCF
    vcf = pysam.VariantFile(ann_vcf_path)
    for rec in vcf:
        ann_info = rec.info.get("ANN")
        if not ann_info:
            continue
        ANN,effects, c_changes, p_changes, proteins, genes = [], [], [], [], [],[]
        
        for ann in ann_info:
            fields = ann.split("|")
            if len(fields) > 1 and fields[1] == "intragenic_variant":
                continue
            ANN.append(ann)
            effects.append(fields[1] if len(fields) > 1 and fields[1] else "-")
            proteins.append(fields[6] if len(fields) > 6 and fields[6] else "-")
            genes.append(fields[3] if len(fields) > 3 and fields[3] else "-")
            c_changes.append(fields[9] if len(fields) > 9 and fields[9] else "-")
            p_changes.append(fields[10] if len(fields) > 10 and fields[10] else "-")
        
        for alt in rec.alts:
            annot_list.append({
                "chrom": rec.chrom,
                "pos": rec.pos,
                "ref": rec.ref,
                "alt": alt,
                "gene": ",".join(genes),
                "protein": ",".join(proteins),
                "effect": ",".join(effects),
                "c_change": ",".join(c_changes),
                "protein_change": ",".join(p_changes),
                "ANN": ",".join(ANN)
            })

columns = ["chrom","pos","ref","alt","gene","protein","effect","c_change","protein_change","ANN"]
annot_df = pd.DataFrame(annot_list, columns=columns,dtype=str)

cov_col = [c for c in ed_df.columns if c.startswith("Coverage-q")][0]
gcov_col = [c for c in ed_df.columns if c.startswith("gCoverage-q")][0]

if "gAllSubs" in ed_df.columns:
    merged_df = annot_df.merge(
        filtered_rows[['Region','Position','AllSubs','Frequency',"sup_reads",'BaseCount[A,C,G,T]',cov_col,'gAllSubs','gFrequency',"gsup_reads",'gBaseCount[A,C,G,T]',gcov_col,"strand"]],
        left_on=['chrom','pos'],
        right_on=['Region','Position'],
        how='left'  
    )
else:
    merged_df = annot_df.merge(
        filtered_rows[['Region','Position','AllSubs','Frequency',"sup_reads",'BaseCount[A,C,G,T]',cov_col,"strand"]],
        left_on=['chrom','pos'],
        right_on=['Region','Position'],
        how='left'  
    )
    
merged_df.drop(columns=['Region','Position'], inplace=True)
#merged_df['BaseCount[A,C,G,T]'] = merged_df['BaseCount[A,C,G,T]'].fillna('-')
#merged_df['Frequency'] = merged_df['Frequency'].fillna('-')
#merged_df['Coverage-q25'] = merged_df['Coverage-q25'].fillna('-')
#merged_df['MeanQ'] = merged_df['MeanQ'].fillna('-')

merged_df["sample"]=sample_name
merged_df["DNA_sample"]=DNA_name
merged_df["spacies"]=snpeff_spacies
merged_df.to_csv(f"{out_dir}{sample_name}_vs_{DNA_name}.edsite_ann.csv",sep="\t",index=False)

