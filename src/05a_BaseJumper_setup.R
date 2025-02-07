#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
bam_dir <- "data/resolveome/DNA/"
fastq_dir <- paste0(wd, "/out/bamtofastq/reads/")

# create bj-wgs samplesheet
# columns: biosampleName, read1, read2
ss <-
  readr::read_csv("out/bamtofastq/samplesheet.csv") %>%
  dplyr::filter(grepl("49686|49882", sample_id)) %>%
  dplyr::transmute(
    biosampleName = sample_id,
    read1 = file.path(fastq_dir, paste0(sample_id, "_1.merged.fastq.gz")),
    read2 = file.path(fastq_dir, paste0(sample_id, "_2.merged.fastq.gz"))) %>%
  dplyr::filter(file.exists(read1), file.exists(read2))
ss %>% readr::write_csv("out/BaseJumper/bj-wgs/samplesheet.csv")

# create bj-somatic-variantcalling samplesheet
# columns: biosampleName, read1, read2, groups, isbulk, bam
ss_tmp <-
  readr::read_tsv("out/nf-resolveome/PD63118/samplesheet.tsv") %>%
  dplyr::filter(run %in% c("49686", "49882")) %>%
  dplyr::transmute(
    biosampleName = id,
    bam = file.path(wd, bam_dir, "WGS", gsub("_.*", "", donor_id),
                    paste0(run, "_", id, ".bam")),
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
ss %>% readr::write_csv("out/BaseJumper/bj-somatic-variantcalling/samplesheet.csv")
