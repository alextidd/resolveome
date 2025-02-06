#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
out_dir <- "out/nf-gatk/"
dir.create(out_dir)

# get samplesheet
readr::read_csv("data/resolveome/DNA/samplesheet_local.csv") %>%
  dplyr::transmute(donor_id, id, bam) %>%
  readr::write_csv(file.path(out_dir, "samplesheet.csv"))
