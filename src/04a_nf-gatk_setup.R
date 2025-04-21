#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
out_dir <- "out/nf-gatk/"
bjsvc_dir <- paste0(wd, "/out/BaseJumper/bj-somatic-variantcalling/dna/")
dir.create(out_dir)

# get sentieon-aligned bams, write samplesheet
readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type == "dna") %>%
  dplyr::transmute(id, donor_id,
                   bam = paste0(bjsvc_dir, donor_id,
                                "_run/secondary_analyses/alignment/",
                                id, ".bam")) %>%
  dplyr::filter(file.exists(bam)) %>%
  readr::write_csv(file.path(out_dir, "samplesheet.csv"))
