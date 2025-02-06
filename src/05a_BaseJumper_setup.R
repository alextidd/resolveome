#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
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
# columns: biosampleName,bam,groups,read1,read2,isbulk
ss <-
  readr::