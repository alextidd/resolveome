#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
fastq_dir <- paste0(wd, "/out/bamtofastq/reads/")
bam_dir <- "data/resolveome/PD63118/"

# read samplesheet
ss <- readr::read_csv("data/resolveome/samplesheet_local.csv")

# create bj-wgs/bj-wes samplesheets
# columns: biosampleName, read1, read2
ss_fastq <-
  readr::read_csv("out/bamtofastq/samplesheet.csv") %>%
  dplyr::transmute(
    biosampleName = sample_id,
    read1 = file.path(fastq_dir, paste0(sample_id, "_1.merged.fastq.gz")),
    read2 = file.path(fastq_dir, paste0(sample_id, "_2.merged.fastq.gz"))) %>%
  dplyr::filter(file.exists(read1), file.exists(read2))
ss_fastq %>%
  dplyr::filter(biosampleName %in% dplyr::filter(ss, seq_type == "dna")$id) %>%
  readr::write_csv("out/BaseJumper/bj-wgs/samplesheet.csv")
ss_fastq %>%
  dplyr::filter(biosampleName %in% dplyr::filter(ss, seq_type == "dnahyb")$id) %>%
  readr::write_csv("out/BaseJumper/bj-wes/samplesheet.csv")

# stage normal bam, create dummy recal table
match_normal <- "/lustre/scratch125/casm/team268im/fa8/117/PTA_49686/PD63118_MATCHED_NORMAL/PD63118.grch38.bam/PD63118.bam"
match_normal_dir <- file.path(wd, "data/lcm/")
match_normal_staged <- paste0(match_normal_dir, "PD63118.bam")
dir.create(match_normal_dir)
system(paste0("ln -s ", match_normal, " ", match_normal_staged))
system(paste0("ln -s ", match_normal, ".bai ", match_normal_staged, ".bai"))
file.create(gsub(".bam", "_recal_data.table", match_normal_staged))

# create bj-somatic-variantcalling samplesheet
ss_pta <-
  ss %>%
  dplyr::filter(seq_type == "dnahyb") %>%
  dplyr::transmute(
    biosampleName = id, groups = donor_id,
    read1 = NA, read2 = NA,
    isbulk = FALSE,
    bam = paste0(wd, "/out/BaseJumper/bj-wes/_250219_235359/secondary_analyses/alignment/", id, ".bam"))
ss_bulk <-
  tibble::tibble(
    id = "PD63118", groups = "PD63118",
    read1 = NA, read2 = NA,
    isbulk = TRUE,
    bam = match_normal_staged
  )
dplyr::bind_rows(ss_pta, ss_bulk) %>%
  dplyr::filter(file.exists(bam)) %>%
  readr::write_csv("out/BaseJumper/bj-somatic-variantcalling/samplesheet.csv")

# create