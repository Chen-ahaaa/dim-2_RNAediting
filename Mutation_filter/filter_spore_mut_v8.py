#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  1 12:41:43 2025

@author: cxn
"""
import pandas as pd
import argparse 
import time
import os
#######################################################################################
start = time.perf_counter() 

#required parameters:--file vcf file  --outpath 
parser = argparse.ArgumentParser()
parser.add_argument("--vcf",type=str,required=True)
parser.add_argument("--outpath",type=str,required=True)

#optional parameters:set min_reads_num 
parser.add_argument("--parents_min_reads",type=int,default=5)
parser.add_argument("--progeny_min_reads",type=int,default=5)
parser.add_argument("--progeny_max_reads",type=int,default=None)
parser.add_argument("--parents_max_reads",type=int,default=None)
parser.add_argument("--case",type=str, nargs="+",default=["case1"])
parser.add_argument("--min_mut_ratio",type=float,default=0.1)
parser.add_argument("--max_parent_mut_reads",type=int,default=5)
parser.add_argument("--max_parent_mut_ratio",type=float,default=0.2)
parser.add_argument("--min_parent_sup_reads",type=int,default=6)
parser.add_argument("--min_progeny_sup_reads",type=int,default=6)
parser.add_argument("--min_sup_ratio",type=float,default=None)
parser.add_argument("--min_sup_strand_reads",type=int,default=1)
parser.add_argument("--only_snp",type=int,default=0)
parser.add_argument("--cross",type=int,default=0)
parser.add_argument("--tetrainfo",type=str,default=None)
parser.add_argument("--dup_file",type=str,default="/mnt/san3/usr/cxn/info/neurospora_genome/newdup/merged_dup_bTH.bed")


args = parser.parse_args()
vcf_file=args.vcf
outpath=args.outpath
parents_min_reads=args.parents_min_reads
progeny_min_reads=args.progeny_min_reads
parents_max_reads=args.parents_max_reads
progeny_max_reads=args.progeny_max_reads
min_parent_sup_reads=args.min_parent_sup_reads
min_progeny_sup_reads=args.min_progeny_sup_reads
min_sup_ratio=args.min_sup_ratio
min_sup_strand_reads=args.min_sup_strand_reads
case=args.case
is_cross=args.cross
tetrainfo=args.tetrainfo
dup_file=args.dup_file
min_mut_ratio=args.min_mut_ratio
max_parent_mut_reads=args.max_parent_mut_reads
max_parent_mut_ratio=args.max_parent_mut_ratio
only_snp=args.only_snp
filename = os.path.basename(vcf_file)
parents_min_reads = parents_min_reads if parents_min_reads is not None else -float("inf")
parents_max_reads = parents_max_reads if parents_max_reads is not None else float("inf")
progeny_min_reads = progeny_min_reads if progeny_min_reads is not None else -float("inf")
progeny_max_reads = progeny_max_reads if progeny_max_reads is not None else float("inf")
#######################################################################################################
#delete the annotation of vcf_file  
def skip_rows_count(vcf_file):
    skip_rows=0
    with open(vcf_file, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('##'):
                skip_rows += 1
            else:
                break
    return skip_rows

def classify_variant(row):
    refs = row["mut_ref"].split(",")
    alts = row["mut"].split(",")
    if all(len(ref) == 1 for ref in refs) and all(len(alt) == 1 for alt in alts):
        return "snp"
    else:
        return "indel"

def add_dp_after_ad(vcf):
    vcf_df=vcf.copy()
    sample_cols = vcf_df.columns[9:]

    first_format = vcf_df.iloc[0]["FORMAT"].split(":")
    
    if "DP" in first_format:
        return vcf_df  
    ad_index = first_format.index("AD") if "AD" in first_format else None

    if ad_index is not None:
        new_format = first_format[:ad_index+1] + ["DP"] + first_format[ad_index+1:]
    else:
        new_format = first_format + ["DP"]
    vcf_df["FORMAT"] = ":".join(new_format)
    
    for col in sample_cols:
        split_vals = vcf_df[col].str.split(":", expand=True)
        
        if ad_index is not None:
            is_dot = split_vals[ad_index] == "."

            ad_nums = split_vals[ad_index].str.split(",", expand=True)
            ad_sum = ad_nums.apply(pd.to_numeric, errors='coerce').sum(axis=1)
            ad_sum = ad_sum.fillna(0).astype(int).astype(str)
            ad_sum[is_dot] = "."

            left = split_vals.iloc[:, :ad_index+1]
            right = split_vals.iloc[:, ad_index+1:]
            split_vals = pd.concat([left, ad_sum.rename("DP"), right], axis=1)
        else:
            split_vals["DP"] = "."
        vcf_df[col] = split_vals.apply(lambda x: ":".join(x), axis=1)
    return vcf_df

def get_mut_RC_SAC(RC_field, mut_idx):
    values = list(map(int, RC_field.split(",")))
    i = mut_idx * 2
    if i+1 < len(values):
        return f"{values[i]},{values[i+1]}"
    else:
        return ".,."


def get_supported_alleles(AD_field, RC_field,
                          min_sup_reads=None, min_sup_ratio=None,
                          min_sup_strand_reads=None):
    if AD_field in (".", None):
        return set(), AD_field

    depths = list(map(int, AD_field.split(",")))
    total = sum(depths)
    if total == 0:
        return set(), AD_field

    alleles = set()
    for i, d in enumerate(depths):
        ratio_ok = (min_sup_ratio is None or (d / total) >= min_sup_ratio)
        reads_ok = (min_sup_reads is None or d >= min_sup_reads)
        if ratio_ok and reads_ok:
            alleles.add(i)

    if min_sup_reads is None and min_sup_ratio is None:
        alleles = {i for i, d in enumerate(depths) if d > 0}

    strand_filtered_out = set()
    if min_sup_strand_reads is not None:
        filtered = set()
        for i in alleles:
            rc_pair = get_mut_RC_SAC(RC_field, i)
            rc1, rc2 = map(int, rc_pair.split(","))
            if rc1 >= min_sup_strand_reads and rc2 >= min_sup_strand_reads:
                filtered.add(i)
        strand_filtered_out = alleles - filtered
        alleles = filtered

    new_depths = []
    for i, d in enumerate(depths):
        if i in strand_filtered_out:
            new_depths.append(".")
        else:
            new_depths.append(str(d))
    new_AD = ",".join(new_depths)

    return alleles, new_AD


def reads_info_by_ad(ad_str,dp_str, mut_idx):
    ad_list = ad_str.split(",")
    if mut_idx >= len(ad_list):
        return ".", "."
    current = ad_list[mut_idx]
    total = int(dp_str)
    if total == 0 or current == ".":
        return ".", "."
    count = int(current)
    ratio = f"{count / total:.2f}"
    return str(count), ratio

def is_het_by_ad(alleles, ad_values, min_ratio=0.3):

    if len(alleles) != 2 or ad_values is None:
        return False
    ad_list = [int(x) for x in ad_values.split(",") if x.isdigit()]
    if len(ad_list) < 2:
        return False
    total = sum(ad_list)
    if total == 0:
        return False
    ratios = [x / total for x in ad_list[:2]]
    return min(ratios) >= min_ratio


def is_hom_by_alleles(alleles, AD_field, min_hom_ratio=0.9):
    if len(alleles) != 1:
        return False
    if AD_field in (".", None):
        return False

    depths = list(map(int, AD_field.split(",")))
    total = sum(depths)
    if total == 0:
        return False

    allele_idx = list(alleles)[0]
    ratio = depths[allele_idx] / total

    return ratio >= min_hom_ratio

def segregation_like_set(p1_alleles, p2_alleles, child_alleles,
                         p1_ad, p2_ad, child_ad, min_ratio=0.3):
    """
    Determine whether a mutation is a "heterozygous segregation type":
      1. One parent is heterozygous (allele frequency ≥ 0.3, 0.3) and the other is homozygous;
         the offspring is homozygous and differs from the homozygous parent;
         the offspring allele corresponds to the other allele of the heterozygous parent.
      2. Both parents are heterozygous (allele frequency ≥ 0.3, 0.3), and the offspring is homozygous.
    """

    p1_is_het = is_het_by_ad(p1_alleles, p1_ad, min_ratio)
    p2_is_het = is_het_by_ad(p2_alleles, p2_ad, min_ratio)
    p1_is_hom = is_hom_by_alleles(p1_alleles,p1_ad)
    p2_is_hom = is_hom_by_alleles(p2_alleles,p2_ad)
    child_is_hom = is_hom_by_alleles(child_alleles,child_ad)

    if not child_is_hom:
        return False

    # one parent heterozygous and the other homozygous
    if (p1_is_het and p2_is_hom) or (p2_is_het and p1_is_hom):
        het_alleles = p1_alleles if p1_is_het else p2_alleles
        hom_alleles = p1_alleles if p1_is_hom else p2_alleles
        if child_alleles != hom_alleles and list(child_alleles)[0] in het_alleles:
            return True

    # both parents heterozygous, offspring homozygous
    if p1_is_het and p2_is_het:
        return True
    return False

def process_site(row):
    tmp = row["REF"] + "," + row["ALT"]
    format_fields = row["FORMAT"].split(":")
    dp_index = format_fields.index("DP") if "DP" in format_fields else None
    ad_index = format_fields.index("AD") if "AD" in format_fields else None
    rc_index = format_fields.index("RC") if "RC" in format_fields else None
    # parents
    parent1_field = row.iloc[9]
    parent2_field = row.iloc[10]

    parent1 = parent1_field.split(":")
    parent2 = parent2_field.split(":")
    DP1 = int(parent1[dp_index]) if len(parent1) > dp_index and parent1[dp_index] != "." else 0
    DP2 = int(parent2[dp_index]) if len(parent2) > dp_index and parent2[dp_index] != "." else 0
    if DP1==0  or DP2==0:
        return []
    if not (parents_min_reads <= DP1 <= parents_max_reads and parents_min_reads <= DP2 <= parents_max_reads):
        return []
    
    alleles1,rAD1 = get_supported_alleles(parent1[ad_index],parent1[rc_index],min_sup_reads=min_parent_sup_reads, min_sup_ratio=min_sup_ratio,min_sup_strand_reads=min_sup_strand_reads)
    alleles2,rAD2 = get_supported_alleles(parent2[ad_index],parent2[rc_index],min_sup_reads=min_parent_sup_reads, min_sup_ratio=min_sup_ratio,min_sup_strand_reads=min_sup_strand_reads)
    if not (len(alleles1)>=1 and len(alleles2)>=1):
        return []
    parent_alleles = alleles1.union(alleles2)
    parent_bases = [tmp.split(",")[a] for a in parent_alleles] if parent_alleles else []
    # progeny
    records = []
    for x in range(11, len(row)):
        sample_field = row.iloc[x]
        
        sample_split = sample_field.split(":")
        DP = int(sample_split[dp_index]) if len(sample_split) > dp_index and sample_split[dp_index] != "." else 0
        if DP==0:
            continue
        if not (progeny_min_reads <= DP <= progeny_max_reads):
            continue

        alleles,rAD = get_supported_alleles(sample_split[ad_index],sample_split[rc_index],min_sup_reads=min_progeny_sup_reads, min_sup_ratio=min_sup_ratio,min_sup_strand_reads=min_sup_strand_reads)

        record_mut = False
        mut_alleles = alleles-parent_alleles 
        mut_case = None
        # have mut
        if len(mut_alleles) > 0:
            mut_ref=",".join(parent_bases) if parent_bases else row["REF"]
            record_mut = True
            if len(mut_alleles) > 1:
                mut_case = "case3"  # multiple mutation directions
            elif len(alleles1) == 1 and len(alleles2) == 1 and len(mut_alleles) == 1:
                mut_case = "case1"  # parents homo and same GT
            else:
                mut_case = "case2"  # other cases: progeny different from parents
        elif  segregation_like_set(alleles1, alleles2, alleles, parent1[ad_index], parent2[ad_index], sample_split[ad_index]):
            record_mut = True
            mut_case = "case5"
            mut_alleles=alleles
            mut_ref=tmp.split(",")[list(parent_alleles - alleles)[0]]
        if record_mut:
            # count mut allele reads & frequency
            for mut_idx in mut_alleles:
                mut_base = tmp.split(",")[mut_idx]
            
                mut_reads_list = []
                for ad, dp in zip([rAD1, rAD2, rAD], [DP1, DP2,  DP]):
                    count, ratio = reads_info_by_ad(ad, dp,mut_idx)
                    mut_reads_list.extend([count, ratio])

                if mut_ref == row["REF"]:
                    mut_origin= "forward"
                elif mut_base == row["REF"]:
                    mut_origin = "reverse"
                else:
                    mut_origin= "complex"
                    
                if len(alleles)==1:
                    homo="homo"
                else:
                    homo=mut_reads_list[-1]
                    
                mut = {
                    "chrom": row["#CHROM"],
                    "pos": row["POS"],
                    "mut_ref": mut_ref,
                    "mut": mut_base,   
                    "mutsample": row.index[x],
                    "origin":mut_origin,
                    "mut_type":"",
                    "homo_type":homo,
                    "mut_reads_ratio": mut_reads_list[-1],
                    "parent_mut_reads": str(mut_reads_list[0])+","+str(mut_reads_list[2]),
                    "mut_case": mut_case
                }
                #if rc_index is not None and len(sample_split) > rc_index and sample_split[rc_index] != ".":
                #    mut_RC = get_mut_RC_SAC(sample_field[rc_index], mut_idx)
                #    mut["mut_FR"] = mut_RC
                
                mut["mut_reads_info"]=",".join(mut_reads_list)
                mut["INFO_AD_DP_RC"]= ":".join([parent1[1],parent1[2],parent1[-1]])+"  "+":".join([parent2[1],parent2[2],parent2[-1]])+"  "+":".join([sample_split[1],sample_split[2],sample_split[-1]])
                records.append(mut)
    return records

def joinsamples(samples):
        return ",".join(samples)

def countsamplenum(samples):
        return len(samples.split(","))

###########################################################################################
annotation_line_num=skip_rows_count(vcf_file)
annotation_lines = ""
with open(vcf_file, 'r') as file:
    for i in range(annotation_line_num):
        line = file.readline().strip()
        annotation_lines += line + '\n'
#import vcf file 
vcf=pd.read_table(vcf_file,sep="\t",skiprows=annotation_line_num)

#filter chrom1~7
vcf=vcf[vcf["#CHROM"].str.contains(r"\.[1-7]$")]

# filter low quality sites
vcf = vcf[~vcf["FILTER"].str.contains("LowQual")].copy()

results = vcf.apply(process_site, axis=1)
flat_records = [rec for site_records in results for rec in site_records]
df = pd.DataFrame(flat_records)

mask = df["parent_mut_reads"].astype(str).str.split(",").apply(
    lambda x: (
        len(x) >= 2
        and all(v.isdigit() and int(v) == 0 for v in x[:2])  # if both are 0
    ) if all(v != "." for v in x[:2]) else False             # if there is a ".", skip
)
df.loc[mask, "homo_type"] = "p" + df.loc[mask, "homo_type"]

#filter sites
df_f=df[df["mut_case"].isin(case)]
df_case5=df_f[df_f["mut_case"]=="case5"]
df_df=df_f[df_f["mut_case"]!="case5"]

mut_reads_split = df_f["mut_reads_info"].astype(str).str.split(",")
p1_mut_freq = mut_reads_split.str[1]
p2_mut_freq = mut_reads_split.str[3]
def valid_freq(x, max_ratio):
    return True if x == "." else float(x) <= max_ratio
mask1 = p1_mut_freq.apply(lambda x: valid_freq(x, max_parent_mut_ratio))
mask2 = p2_mut_freq.apply(lambda x: valid_freq(x, max_parent_mut_ratio))
df_f = df_f[mask1 & mask2]

df_f = df_f[
    df_f["parent_mut_reads"].apply(
        lambda x: all(
            v == "." or (v.isdigit() and int(v) <= max_parent_mut_reads)
            for v in x.split(",")))]
df_f = pd.concat([df_f, df_case5], ignore_index=True)

##
df_f["mut_type"]=df_f.apply( classify_variant,axis=1)
df_f["mut_direction"]=df_f["mut_ref"]+">"+df_f["mut"]
df_f.loc[df_f["mut_type"] == "snp", "mut_direction"] = df_f.loc[df_f["mut_type"] == "snp", "mut_direction"].replace({"G>A": "C>T", "A>G": "T>C", "A>T": "T>A", "G>T": "C>A", "A>C": "T>G", "G>C": "C>G"})

if only_snp==1:
    df_f=df_f[df_f["mut_type"]=="snp"]

df_f["spore_num"]=len(vcf.columns)-11
group=filename.split(".")[0]
df_f["group"]=group

comment_lines = [
    f"# [vcf_file]: {vcf_file}",
    "# [chrom]; [pos]; [mut_ref]: parents allele; [mut]: mutant sample allele; [mutsample]",
    "# [origin]: forward--mutant site is the same as reference before the mutation; reverse--mutant site is the same as reference after the mutation; complex--other situations",
    "# [mut_type]: snp or indel",
    "# [homo_type]: phomo--parents have no reads supporting mutant allele and progeny is homozygous mutation; homo--parents have reads supporting mutant allele and progeny is homozygous mutation; otherwise, mark the mut_reads_ratio)",
    "# [mut_reads_ratio]: the proportion of reads supporting the mutant allele in the mutant sample",
    "# [parent_mut_reads]: the number of reads supporting the mutant allele in parents,dot indicates the supporting reads > min_sup but in a single direction",
    "# [mut_case]: case1--parents are homozygous and have the same GT, progeny has one type different allele(not appears in parents); case2--other possible parental conditions in case1; case3--multiple-directional mutations; case4--the parents are homozygous and have different GT, all the progenies(within one cross) have the same GT as one of the parents; case5--heterozygous parents",
    "# [mut_reads_info]: the number and proportion of reads supporting the mutant allele in the parents and progeny (parent1_count, parent1_ratio; parent2_count, parent2_ratio; progeny_count, progeny_ratio),dot indicates the supporting reads > min_sup but in a single direction",
    "# [INFO_AD_DP_RC]: parent1  parent2  mutsample",
    "# [mut_direction]: snp--6 types, indel--mut_ref>mut",
    "# [spore_num]: number of all sequenced spores within the [group]",
    f"# *filter standard*: [case]--{case}; [parents_min_reads]--{parents_min_reads}; [parents_max_reads]--{parents_max_reads}; [max_parent_mut_ratio]--{max_parent_mut_ratio}; [progeny_min_reads]--{progeny_min_reads}; [progeny_max_reads]--{progeny_max_reads}; [min_mut_ratio]--{min_mut_ratio}, [max_parent_mut_reads]--{max_parent_mut_reads}; [only_snp]--{only_snp}; [min_parent_sup_reads]--{min_parent_sup_reads}; [min_progeny_sup_reads]--{min_progeny_sup_reads}; [min_sup_ratio]--{min_sup_ratio}; [min_sup_strand_reads]--{min_sup_strand_reads}"
]
# "# [mut_FR]: the number of reads supporting mutant allele on the forward and reverse strands of mutsample",

filter_file = f"{outpath}{filename}.filter.csv"

with open(filter_file, "w") as f:
    for line in comment_lines:
        f.write(line + "\n")
    df_f.to_csv(f,sep="\t",index=False)

##################################################################################
#####################################################################################
def filter_case4(row,cross_dict):
    tmp = row["REF"] + "," + row["ALT"]
    format_fields = row["FORMAT"].split(":")
    dp_index = format_fields.index("DP") if "DP" in format_fields else None
    ad_index = format_fields.index("AD") if "AD" in format_fields else None
    rc_index = format_fields.index("RC") if "RC" in format_fields else None
    # parents
    parent1_field = row.iloc[9]
    parent2_field = row.iloc[10]
    parent1 = parent1_field.split(":")
    parent2 = parent2_field.split(":")
    DP1 = int(parent1[dp_index]) if len(parent1) > dp_index and parent1[dp_index] != "." else 0
    DP2 = int(parent2[dp_index]) if len(parent2) > dp_index and parent2[dp_index] != "." else 0
    if DP1==0  or DP2==0:
        return []
    if not (parents_min_reads <= DP1 <= parents_max_reads and parents_min_reads <= DP2 <= parents_max_reads):
        return []
    alleles1,rAD1= get_supported_alleles(parent1[ad_index],parent1[rc_index],min_sup_reads=min_parent_sup_reads, min_sup_ratio=min_sup_ratio,min_sup_strand_reads=min_sup_strand_reads)
    alleles2,rAD2 = get_supported_alleles(parent2[ad_index],parent2[rc_index],min_sup_reads=min_parent_sup_reads, min_sup_ratio=min_sup_ratio,min_sup_strand_reads=min_sup_strand_reads)

    records = []
    # parent homo and same GT
    if (len(alleles1) == 1) and (len(alleles2) == 1) and (alleles1 != alleles2):
        for cross, spores in cross_dict.items():
            spore_gts = []
            valid_group = True
            for spore in spores:
                if spore not in vcf.columns:
                    continue
                sample_field=row[spore]
                sample_split= sample_field.split(":")
                DP=int(sample_split[dp_index]) if len(sample_split) > dp_index and sample_split[dp_index] != "." else 0
                if DP==0:
                    valid_group = False
                    break
                alleles,rAD = get_supported_alleles(sample_split[ad_index],sample_split[rc_index],min_sup_reads=min_progeny_sup_reads, min_sup_ratio=min_sup_ratio,min_sup_strand_reads=min_sup_strand_reads)
                if not ((progeny_min_reads <= DP <= progeny_max_reads) and len(alleles)==1):
                    valid_group = False
                    break
                spore_gts.extend(list(alleles))
            if not valid_group or not spore_gts:
                continue
            if len(set(spore_gts)) == 1:
                if set(spore_gts)==alleles1 :
                    mut_ref =  tmp.split(",")[next(iter(alleles2))]
                    mut_base= tmp.split(",")[next(iter(alleles1))]
                    cross_mut_parent = "B"
                elif set(spore_gts)==alleles2 :
                    mut_ref =  tmp.split(",")[next(iter(alleles1))]
                    mut_base= tmp.split(",")[next(iter(alleles2))]
                    cross_mut_parent = "A"
                    
                if mut_ref == row["REF"]:
                    mut_origin= "forward"
                elif mut_base == row["REF"]:
                    mut_origin = "reverse"
                else:
                    mut_origin= "complex"
                    
                records.append({
                        "chrom": row["#CHROM"],
                        "pos": row["POS"],
                        "mut_ref": mut_ref,
                        "mut": mut_base,
                        "origin":mut_origin,
                        "mut_cross": cross,
                        "cross_mut_samples": "uncertain",
                        "mut_case": "case4",
                        "cross_mut_parent": cross_mut_parent,
                        "mut_tetra_type": "2:2",
                        "tetra_parent_type": "2:2",
                        "cross_mut_ratio": "homo",
                        "parent_mut_ratio": "homo",
                        "homo_type": "homo",
                        "cross_mut_INFO": "",
                        "sequenced_spore":len(spores),
                        "spore_num":len(vcf.columns)-11,
                        "cross_num":len(cross_dict)
                    })
    return records
    

############################################################################
if is_cross==1 and tetrainfo is not None:
    tetra=pd.read_csv(tetrainfo,comment="#",sep="\t",header=None)
    cross_dict = {}
    for _, row in tetra.iterrows():
        cross = row.iloc[0]
        spores = []
        for col in row.iloc[1:]:
            if pd.isna(col):
                continue
            spores.extend(col.split(","))
        cross_dict[cross] = spores
    #add 2:2>4:0 mut
    case4_results = vcf.apply(filter_case4, cross_dict=cross_dict, axis=1)
    flat_records = [rec for site_records in case4_results for rec in site_records]
    case4 = pd.DataFrame(flat_records)
    case4["mut_type"]=case4.apply(classify_variant,axis=1)
    case4["mut_direction"]=case4["mut_ref"]+">"+case4["mut"]
    case4.loc[case4["mut_type"] == "snp", "mut_direction"] = case4.loc[case4["mut_type"] == "snp", "mut_direction"].replace({"G>A": "C>T", "A>G": "T>C", "A>T": "T>A", "G>T": "C>A", "A>C": "T>G", "G>C": "C>G"})
    case4["group"]=group
    if only_snp==1:
        case4=case4[case4["mut_type"]=="snp"]
    case4_file = f"{outpath}{filename}.case4.csv"
    case4.to_csv(case4_file,sep="\t",index=False)
        
end = time.perf_counter()
print(f"finished {vcf_file} [time: {end - start:.2f}s]")















