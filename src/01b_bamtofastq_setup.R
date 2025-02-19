#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
dir.create("out/bamtofastq/")

# generate samplesheet
readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::transmute(mapped = bam, index = paste0(bam, ".bai"),
                   sample_id = id, file_type = "bam") %>%
  readr::write_csv("out/bamtofastq/samplesheet.csv")
