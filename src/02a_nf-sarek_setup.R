#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
out_dir <- "out/nf-sarek/"
dir.create(out_dir)

# get samplesheet
ss <- readr::read_csv("data/resolveome/DNA/samplesheet_local.csv")

# generate sarek samplesheet
ss %>%
  dplyr::transmute(patient = donor_id, sample = id, bam,
                   bai = paste0(bam, ".bai")) %>%
  readr::write_csv(file.path(out_dir, "samplesheet.csv"))
