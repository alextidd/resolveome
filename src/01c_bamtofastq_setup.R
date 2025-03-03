#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
dir.create("out/bamtofastq/", showWarnings = FALSE)

# get samplesheet
ss <- readr::read_csv("data/resolveome/samplesheet_local.csv")

# generate bamtofastq samplesheet
ss_unmerged <-
  ss %>% dplyr::filter(!run_id %in% c(49901, 50072)) %>%
  dplyr::transmute(sample_id = id, mapped = bam, index = paste0(bam, ".bai"),
                   file_type = "bam")

# generate bamtofastq samplesheet - must revise bam paths for merged samples
ss_merged <-
  ss %>% dplyr::filter(run_id %in% c(49901, 50072)) %>%
  dplyr::mutate(sample_id = paste0(cell_id, "_", seq_type, "_merged"),
                mapped = paste0(wd, "/out/merge_bams/", donor_id, "/",
                                sample_id, "/bam/", sample_id, ".bam"),
                index = paste0(mapped, ".bai"),
                file_type = "bam") %>%
  dplyr::distinct(sample_id, mapped, index, file_type)

# write samplesheet
ss_merged %>%
  write_csv("out/bamtofastq/samplesheet_merged.csv")
dplyr::bind_rows(ss_unmerged, ss_merged) %>%
  readr::write_csv("out/bamtofastq/samplesheet.csv")
