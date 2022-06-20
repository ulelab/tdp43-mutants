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

def print_comparisons(snakemake_files, nf_quantseq_files):
    for snakemake_file_name in snakemake_files:
        if snakemake_file_name.endswith(".gz"):
            in_f = gzip.open(snakemake_file_name, "rt")
        else:
            in_f = open(snakemake_file_name)
        try:
            snakemake_hash = hashlib.md5(in_f.read().encode()).hexdigest()
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
            nf_quantseq_hash = hashlib.md5(in_f.read().encode()).hexdigest()
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

# Cutadapt
print("--- {} ---".format("cutadapt"))
print_comparisons(
    glob.glob(snakemake_results + "cutadapt/*"),
    glob.glob(nf_quantseq_results + "cutadapt/*")
)

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

# Remove random priming

print("--- {} ---".format("Remove random priming"))
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.bed"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.bed")
)
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.unique.bed"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.unique.bed")
)
print_comparisons(
    glob.glob(snakemake_results + "mergedclusters/polyAclusters.unique.annotated.bed"),
    glob.glob(nf_quantseq_results + "polyaclusters/merged_polya.unique.annotated.bed")
)

# Count tables

print("--- {} ---".format("Count tables"))
print_comparisons(
    glob.glob(snakemake_results + "counts/*bed.gz"),
    glob.glob(nf_quantseq_results + "bedtools/*window_merged_polya.bed")
)
