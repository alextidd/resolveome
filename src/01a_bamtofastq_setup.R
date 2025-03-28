#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
dir.create("out/bamtofastq/", showWarnings = FALSE)

# write samplesheet
readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::transmute(sample_id = id, mapped = bam, index = paste0(bam, ".bai"),
                   file_type = "bam") %>%
  readr::write_csv("out/bamtofastq/samplesheet.csv")
