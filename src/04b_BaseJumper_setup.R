#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
out_dir <- "out/BaseJumper/"

# create samplesheet
# columns: biosampleName, read1, read2, groups, isbulk, bam
ss <-
  readr::read_tsv("out/nf-resolveome/PD63118/samplesheet.tsv") %>%
  dplyr::transmute(
    biosampleName = id,
    bam = file.path(wd, out_dir, donor_id, id, "bam", paste0(id, ".cram")),
    groups = donor_id, read1 = NA, read2 = NA,
    isbulk = FALSE)

# create bulk entries
ss_bulk <-
  ss %>%
  dplyr::distinct(groups, read1, read2) %>%
  dplyr::mutate(
    biosampleName = "normal",
    bam = "/nfs/irods-cgp-sb10-sdb/intproj/3464/sample/PD63118b_lo0001/PD63118b_lo0001.v1.sample.dupmarked.bam",
    isbulk = TRUE)

# combine
ss <- dplyr::bind_rows(ss, ss_bulk)
ss %>% readr::write_csv("out/BaseJumper/samplesheet.csv")