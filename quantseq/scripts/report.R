#!/usr/bin/env Rscript

pwd = getwd()
if(!dir.exists("results/plots")) dir.create("results/plots")

params = list(pas = snakemake@input[["pas"]], polyasite2 = snakemake@params[["polyasite2"]], rds = snakemake@input[["rds"]])

print(params)

rmarkdown::render("scripts/report.Rmd", 
    params = params, 
    output_file = snakemake@output[["report"]],
    output_dir = file.path(pwd, "results"),
    knit_root_dir = pwd)