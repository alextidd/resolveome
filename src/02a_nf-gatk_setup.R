#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
out_dir <- "out/nf-gatk/"
test_dir <- "out/test/nf-gatk/"
dir.create(out_dir)
dir.create(test_dir)

# write samplesheet
readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type == "dna") %>%
  dplyr::transmute(donor_id, id, bam) %>%
  readr::write_csv(file.path(out_dir, "samplesheet.csv"))

# write test samplesheet
readr::read_csv("out/test/samplesheet.csv") %>%
  dplyr::filter(seq_type == "dna") %>%
  dplyr::transmute(donor_id, id, bam) %>%
  readr::write_csv("out/test/nf-gatk/samplesheet.csv")
