# libraries
library(magrittr)

# dirs
dir.create("out/vdj_coverage/regions/", showWarnings = FALSE)
dir.create("out/vdj_coverage/all_cells/", showWarnings = FALSE)
dir.create("out/vdj_coverage/filter_cells/", showWarnings = FALSE)

# load data
df <-
  readr::read_tsv("data/vdj_coverage/ig_tcr_genes_pseudogenes.tsv") %>%
  janitor::clean_names() %>%
  dplyr::select(chr = chromosome_scaffold_name, start = gene_start_bp,
                end = gene_end_bp, gene = gene_name, everything()) %>%
  dplyr::filter(chr %in% c("2", "7", "14", "22")) %>%
  dplyr::mutate(class = dplyr::case_when(grepl("^IG", gene) ~ "BCR",
                                         grepl("^TR", gene) ~ "TCR",
                                         TRUE ~ NA))

# write tsv of genes
df %>%
  readr::write_tsv("out/vdj_coverage/regions/ig_tcr_genes.tsv")

# write bed of genes
df %>%
  dplyr::select(chr, start, end, gene) %>%
  readr::write_tsv("out/vdj_coverage/regions/ig_tcr_genes.bed",
                   col_names = FALSE)

# summarise BCR/TCR regions per chr
df_regions <-
  df %>%
  dplyr::group_by(chr, class) %>%
  dplyr::summarise(start = min(start), end = max(end)) %>%
  dplyr::select(chr, start, end, class)

# write bed of regions
df_regions %>%
  dplyr::select(chr, start, end, class) %>%
  readr::write_tsv("out/vdj_coverage/regions/ig_tcr_regions.bed",
                   col_names = FALSE)

# get samplesheet of all cells
readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type == "dna") %>%
  readr::write_csv("out/vdj_coverage/all_cells/samplesheet.csv")

# get cell filters from manual inspection
man_insp <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::filter((!chr_dropout & !suspected_doublet) | plate == 10)

# write samplesheet of filtered cells
readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type %in% c("dna", "dnahyb"),
                cell_id %in% man_insp$cell_id) %>%
  readr::write_csv("out/vdj_coverage/filter_cells/samplesheet.csv")
