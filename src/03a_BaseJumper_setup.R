#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
fastq_dir <- paste0(wd, "/out/bamtofastq/reads/")

# read samplesheet
ss <- readr::read_csv("data/resolveome/samplesheet_local.csv")

# read manual inspection results, get clean cell ids (include all of plate 10)
clean_cell_ids <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::filter((!suspected_doublet & !chr_dropout) | plate == 10) %>%
  dplyr::pull(cell_id) %>%
  unique()

# create dnahyb/dna bj-somatic-variantcalling samplesheets
# columns: biosampleName,read1,read2,groups,isbulk,bam
ss_fastq <-
  ss %>%
  dplyr::mutate(filter_lvl = "all_cells") %>%
  dplyr::bind_rows(ss %>% dplyr::filter(cell_id %in% clean_cell_ids) %>% dplyr::mutate(filter_lvl = "filter_cells")) %>%
  dplyr::transmute(
    biosampleName = id,
    read1 = file.path(fastq_dir, paste0(id, "_1.merged.fastq.gz")),
    read2 = file.path(fastq_dir, paste0(id, "_2.merged.fastq.gz")),
    groups = donor_id, isbulk = FALSE, bam = "", seq_type, filter_lvl) %>%
  {split(., .$filter_lvl)} %>%
  purrr::map(~ split(.x %>% dplyr::select(-seq_type, -filter_lvl), .x$seq_type))

# save all cells for bj-dna-qc dna
ss_fastq$all_cells$dna %>%
  dplyr::select(biosampleName, read1, read2) %>%
  readr::write_csv("out/BaseJumper/bj-dna-qc/dna/samplesheet.csv")

# save filtered cells for bj-somatic-variantcalling dna
ss_fastq$filter_cells$dna %>%
  readr::write_csv("out/BaseJumper/bj-somatic-variantcalling/dna/samplesheet.csv")

# save filtered cells for bj-somatic-variantcalling dnahyb
ss_fastq$filter_cells$dnahyb %>%
  readr::write_csv("out/BaseJumper/bj-somatic-variantcalling/dnahyb/samplesheet.csv")

# save all cells for bj-expression rna
ss_fastq$all_cells$rna %>%
  dplyr::select(biosampleName, read1, read2) %>%
  readr::write_csv("out/BaseJumper/bj-expression/rna/samplesheet.csv")