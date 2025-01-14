#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
bam_dir <- "data/resolveome/"

# create samplesheet
# columns: biosampleName, read1, read2, groups, isbulk, bam
ss_tmp <-
  readr::read_tsv("out/nf-resolveome/PD63118/samplesheet.tsv") %>%
  dplyr::transmute(
    biosampleName = id,
    bam = file.path(wd, bam_dir, donor_id, id, "bam", paste0(id, ".bam")),
    groups = donor_id, read1 = NA, read2 = NA,
    isbulk = FALSE)

# create bulk entries
ss_bulk <-
  ss_tmp %>%
  dplyr::distinct(groups, read1, read2) %>%
  dplyr::mutate(
    biosampleName = "normal",
    bam = "/nfs/irods-cgp-sb10-sdb/intproj/3464/sample/PD63118b_lo0001/PD63118b_lo0001.v1.sample.dupmarked.bam",
    isbulk = TRUE)

# combine
ss <- dplyr::bind_rows(ss_tmp, ss_bulk)
ss %>% readr::write_csv("out/BaseJumper/samplesheet.csv")
