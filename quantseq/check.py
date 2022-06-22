#!/usr/local/bin/python
# -*- coding: utf-8 -*-

import glob
import hashlib
import gzip
import os
import difflib
import sys

###############
### GLOBALS ###
###############

char_count_to_show = 40

#################
### FUNCTIONS ###
#################

def compare_str(str1, str2):
    s = difflib.SequenceMatcher(None, str1, str2)
    return(s.find_longest_match(0, len(str1), 0, len(str2)).size)

def get_similar_file_name(str_in, list_in):
    scores = [
        (
            compare_str(
                os.path.basename(str_in),
                os.path.basename(str_comparison)
            ),
            str_comparison
        )
        for str_comparison in list_in
    ]
    return(sorted(scores)[-1][1])

def print_comparisons(snakemake_files, nf_quantseq_files, presort_lines=False):
    for snakemake_file_name in snakemake_files:
        if snakemake_file_name.endswith(".gz"):
            in_f = gzip.open(snakemake_file_name, "rt")
        else:
            in_f = open(snakemake_file_name)
        try:
            in_f_str = in_f.read()
            if presort_lines:
                in_f_str = "\n".join(sorted(in_f_str.split("\n")))
            snakemake_hash = hashlib.md5(in_f_str.encode()).hexdigest()
        except UnicodeDecodeError:
            sys.exit("Unicode error reading {}".format(snakemake_file_name))
        in_f.close()
        nf_quantseq_file_name = get_similar_file_name(
            snakemake_file_name,
            nf_quantseq_files
        )
        if nf_quantseq_file_name.endswith(".gz"):
            in_f = gzip.open(nf_quantseq_file_name, "rt")
        else:
            in_f = open(nf_quantseq_file_name)
        try:
            in_f_str = in_f.read()
            if presort_lines:
                in_f_str = "\n".join(sorted(in_f_str.split("\n")))
            nf_quantseq_hash = hashlib.md5(in_f_str.encode()).hexdigest()
        except UnicodeDecodeError:
            sys.exit("Unicode error reading {}".format(nf_quantseq_file_name))
        in_f.close()
        print("{} - Comparing {} ({}) and {} ({})".format(
            "✅" if snakemake_hash == nf_quantseq_hash else "❌",
            "..." + snakemake_file_name[
                len(snakemake_file_name)-char_count_to_show:len(snakemake_file_name)
            ],
            snakemake_hash,
            "..." + nf_quantseq_file_name[
                len(nf_quantseq_file_name)-char_count_to_show:len(nf_quantseq_file_name)
            ],
            nf_quantseq_hash,
        ))

############
### MAIN ###
############

snakemake_results = "/camp/home/jonesm5/home/users/jonesm5/projects/tdp43-mutants/quantseq/results/"
nf_quantseq_results = "/camp/home/jonesm5/home/users/jonesm5/projects/nf-quantseq/results/"

### cutadapt directory
print("=== {} ===".format("cutadapt"))

# Cutadapt
print("--- {} ---".format("cutadapt"))
print_comparisons(
    glob.glob(snakemake_results + "cutadapt/*"),
    glob.glob(nf_quantseq_results + "cutadapt/*")
)

### polya directory
print("=== {} ===".format("polya"))

# Cutadapt PolyA trimmed
print("--- {} ---".format("cutadapt polyA trimmed"))
print_comparisons(
    glob.glob(snakemake_results + "polya/*trimmed*"),
    glob.glob(nf_quantseq_results + "polya_cutadapt/*fastq*")
)

# PolyA Clusters Merged

print("--- {} ---".format("polyA clusters merged"))
print_comparisons(
    glob.glob(snakemake_results + "polya_mapped/merged/*bedgraph"),
    glob.glob(nf_quantseq_results + "polyaclusters/*bedgraph")
)

### mergedclusters directory
print("=== {} ===".format("mergedclusters"))

# Remove random priming

print("--- {} ---".format("Remove random priming"))
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.bed"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.bed")
)
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.bedgraph"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.bedgraph")
)
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.unique.bed"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.unique.bed")
)
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.unique.bedgraph"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.unique.bedgraph")
)
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.unique.annotated.bed"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.unique.annotated.bed")
)

### mapped directory
print("=== {} ===".format("mapped"))

print_comparisons(
    glob.glob(snakemake_results + "mapped/*fastq.gz"),
    glob.glob(nf_quantseq_results + "mapped/*fastq.gz")
)

### counts directory
print("=== {} ===".format("counts"))

# bamtobed

print("--- {} ---".format("bamtobed"))
print_comparisons(
    glob.glob(snakemake_results + "counts/*Aligned.sortedByCoord.out.bed"),
    glob.glob(nf_quantseq_results + "bedtools/*polya_trimmed.bed"),
    presort_lines=True
)

# Count tables

print("--- {} ---".format("Count tables"))
print_comparisons(
    glob.glob(snakemake_results + "counts/*bed.gz"),
    glob.glob(nf_quantseq_results + "bedtools/*window_merged_polya.bed")
)

# Bedgraphs

print("--- {} ---".format("bedgraph"))
print_comparisons(
    glob.glob(snakemake_results + "counts/*bedgraph.gz"),
    glob.glob(nf_quantseq_results + "bedgraphconvert/*")
)
