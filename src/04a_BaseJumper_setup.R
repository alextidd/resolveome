#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
fastq_dir <- paste0(wd, "/out/bamtofastq/reads/")

# read samplesheet
ss <- readr::read_csv("data/resolveome/samplesheet_local.csv")

# create wes/wgs bj-somatic-variantcalling samplesheets
# columns: biosampleName,read1,read2,groups,isbulk,bam
ss_fastq <-
  ss %>%
  dplyr::transmute(
    biosampleName = id,
    read1 = file.path(fastq_dir, paste0(id, "_1.merged.fastq.gz")),
    read2 = file.path(fastq_dir, paste0(id, "_2.merged.fastq.gz")),
    groups = donor_id, isbulk = FALSE, bam = "",
    seq_type) %>%
  dplyr::filter(file.exists(read1), file.exists(read2)) %>%
  {split(., .$seq_type)} %>%
  lapply(dplyr::select, -seq_type)

# save for bj-dna-qc wgs
ss_fastq$dna %>%
  dplyr::select(biosampleName, read1, read2) %>%
  readr::write_csv("out/BaseJumper/bj-dna-qc/wgs/samplesheet.csv")

# save for bj-somatic-variantcalling wgs
ss_fastq$dna %>%
  readr::write_csv("out/BaseJumper/bj-somatic-variantcalling/wgs/samplesheet.csv")

# save for bj-somatic-variantcalling wes
ss_fastq$dnahyb %>%
  readr::write_csv("out/BaseJumper/bj-somatic-variantcalling/wes/samplesheet.csv")

# create bj-expression samplesheet for rna
ss_fastq$rna %>%
  dplyr::select(biosampleName, read1, read2) %>%
  readr::write_csv("out/BaseJumper/bj-expression/samplesheet.csv")
