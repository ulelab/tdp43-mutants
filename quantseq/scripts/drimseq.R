#!/usr/bin/env Rscript

# Script to run DRIMseq analysis with batch modelling
# A. M. Chakrabarti
# Last updated: 15th January 2021

library(data.table)
library(DRIMSeq)
library(BiocParallel)
library(tictoc)

args <- commandArgs(trailingOnly = TRUE)
print(args)

mutant <- args[1]
pas <- args[2]
counttables <- args[3]
counttables <- strsplit(counttables, ",")[[1]]
output <- args[4]
threads <- as.integer(args[5])

message(mutant)

# ==========
# Load PAS data
# ==========
message("Loading PAS data")

pas.dt <- fread(pas, col.names = c("seqnames", "start", "end", "id", "score", "strand", "ensg", "hgnc", "region"))
pas.dt <- pas.dt[region == "UTR3"]
pas.dt[, total_score := sum(score), by = ensg]
pas.dt[, percent_5 := 0.05 * total_score]
pas.dt[, percent_1 := 0.01 * total_score]
pas.dt[, multi := .N, by = ensg]

# ==========
# Load count data
# ==========
message("Loading count data")

# Get file paths and remove dud experiments
# all.count.files <- list.files("results/counts", full.names = TRUE, pattern = ".bed.gz$")
all.count.files <- counttables
count.files <- grep("WT_neg_siTDP1_2|M337P_pos_siTDP1_2|G335A_neg_siTDP1_1|G294A_pos_siTDP1_1|A326P_neg_siTDP1_3", all.count.files, invert = TRUE, value = TRUE)

# Load and merge all count files
count.list <- lapply(1:length(count.files), function(i) fread(count.files[[i]], select = c(4, 7, 8, 10), col.names = c("PAS", "ensg", "hgnc", gsub(".polyacount.bed.gz", "", basename(count.files)[i]))))
count.dt <- Reduce(function(x, y) merge(x = x, y = y, by = c("PAS", "ensg", "hgnc")), count.list)
rm(count.list)

# Remove intergenic PAS and those genes with only 1 PAS
count.dt <- count.dt[ensg != "intergenic"]
count.dt <- count.dt[PAS %in% pas.dt[multi > 1]$id]

# ==========
# Run DRIMseq
# ==========
message("Running DRIMseq")

tic()

BPPARAM = MulticoreParam(workers=threads)

counts.df <- as.data.frame(count.dt[, `:=` (gene_id = paste0(hgnc, "_", ensg),
                                            feature_id = PAS)])
counts.df <- counts.df[, grep(paste0(mutant, "|gene_id|feature_id"), colnames(counts.df))]
if(mutant == "WT") counts.df <- counts.df[, grep("WTssIV", colnames(counts.df), invert = TRUE)]

samples.df <- data.frame(sample_id = head(colnames(counts.df), -2))
samples.df$group <- factor(ifelse(grepl("pos", samples.df$sample_id), "Rescue", "KD"), levels = c("KD", "Rescue"))
samples.df$mutant <- sapply(strsplit(as.character(samples.df$sample_id), "\\_"), "[[", 1)
samples.df$batch <- sapply(strsplit(as.character(samples.df$sample_id), "\\_"), "[[", 3)

print(samples.df$sample_id)
message(floor(nrow(samples.df) * 0.75))
message(floor(min(table(samples.df$group)) * 0.75))

d <- dmDSdata(counts = counts.df, samples = samples.df)

# Require gene expression in at least 75% of all samples, and PAS expression in 0.75% of samples for either KD or rescue (whichever has fewest)
d <- dmFilter(d, 
            min_samps_gene_expr = floor(nrow(samples.df) * 0.75), 
            min_samps_feature_expr = floor(min(table(samples.df$group)) * 0.75), 
            min_gene_expr = 10, 
            min_feature_expr = 5)

if(mutant %in% c("G294A", "Q331K", "WTssIV")) {
design_full <- model.matrix(~ group, data = DRIMSeq::samples(d))
} else {
design_full <- model.matrix(~ group + batch, data = DRIMSeq::samples(d))
}

set.seed(42)
d <- dmPrecision(d, design = design_full, verbose = 1, BPPARAM = BPPARAM)
d <- dmFit(d, design = design_full, verbose = 1, BPPARAM = BPPARAM)
d <- dmTest(d, coef = "groupRescue", verbose = 1, BPPARAM = BPPARAM)

# if(!dir.exists("drimseq")) dir.create("drimseq")

# saveRDS(d, file = paste0("drimseq/", mutant, ".d.rds"))
saveRDS(d, file = output)

toc()