# libraries
library(magrittr)
library(ggplot2)
library(patchwork)

# dirs


# load samplesheet
ss <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::mutate(run_id = as.character(run_id))
ss_exp <-
  ss %>%
  dplyr::select(id, cell_id, seq_type) %>%
  dplyr::group_by(id, cell_id) %>%
  tidyr::expand_grid(name = c("celltype_SHM", "celltype_VDJ_recomb", "celltype_bj")) %>%
  dplyr::ungroup()

# load annots
annots <-
  "data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv" %>%
  readr::read_tsv()

# load read counts
read_counts <-
  "out/BaseJumper/bj-expression/PD63118_250225_102304/read_counts/combined_read_counts.txt" %>%
  readr::read_csv() %>%
  dplyr::rename(id = sample_name) %>%
  dplyr::left_join(ss) %>%

# plot read counts
read_counts %>%
  ggplot(aes(x = reorder(id, -reads), y = reads, fill = run_id)) +
  geom_col() +
  guides(x = guide_axis(angle = -90)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  labs(x = "")

# compare the same cells between the two runs
read_counts %>%
  dplyr::group_by(cell_id) %>%
  dplyr::mutate(max_reads = max(reads)) %>%
  dplyr::filter(plate == 3) %>%
  ggplot(aes(x = reorder(cell_id, -max_reads), y = reads, fill = run_id)) +
  geom_col(position = "dodge") +
  guides(x = guide_axis(angle = -90)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "")

# get celltype curated markers from celltypist
ct_markers <-
  readr::read_tsv("../reference/celltypist/Basic_celltype_information.tsv") %>%
  janitor::clean_names() %>%
  tidyr::separate_longer_delim("curated_markers", delim = ", ")

# load gene counts
gene_counts <-
  "out/BaseJumper/bj-expression/PD63118_250225_102304/secondary_analyses/quantification_salmon/df_gene_counts_salmon.tsv" %>%
  readr::read_tsv() %>% 
  dplyr::filter(gene_symbol %in% ct_markers$curated_markers) %>%
  dplyr::mutate(id = gsub("salmon_outdir_", "", File)) %>%
  dplyr::left_join(ss)

# get n genes and n reads per id
n_genes_per_id <-
  gene_counts %>%
  dplyr::filter(countsFromAbundanceNo > 0) %>%
  dplyr::group_by(id, run_id) %>%
  dplyr::summarise(n = dplyr::n_distinct(gene_symbol_gene_id), name = "n genes")
n_reads_per_id <-
  read_counts %>%
  dplyr::transmute(id, run_id, n = reads, name = "n reads")

p_dat <-
  dplyr::bind_rows(n_genes_per_id, n_reads_per_id) %>%
  dplyr::left_join(n_reads_per_id %>% dplyr::transmute(id, n_reads = n)) %>%
  dplyr::left_join(ss)

# plot all runs
p_dat %>%
  ggplot(aes(x = reorder(id, -n_reads), y = n, fill = run_id)) +
  geom_col() +
  guides(x = guide_axis(angle = -90)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "") +
  facet_wrap(~ name, scales = "free_y", ncol = 1) +
  theme_classic()

# compare resequencing
p_dat %>%
  dplyr::filter(plate == 3) %>%
  ggplot(aes(x = cell_id, y = n, fill = run_id)) +
  geom_col(position = "dodge") +
  guides(x = guide_axis(angle = -90)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "") +
  facet_wrap(~ name, scales = "free_y", ncol = 1) +
  theme_classic()

# convert to matrix for heatmap
gene_counts %>%
  dplyr::select(id, gene_symbol_gene_id,
                count = countsFromAbundancelengthScaledTPM) %>%
  tidyr::pivot_wider(names_from = gene_symbol_gene_id, values_from = count) %>%
  tibble::column_to_rownames("id") %>%
  as.matrix()

# load celltype predictions from inspection of VDJ recombination
cts_vdj <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::select(id, celltype_SHM, celltype_VDJ_recomb, suspected_doublet) %>%
  tidyr::pivot_longer(cols = c("celltype_SHM", "celltype_VDJ_recomb"))

# load basejumper celltype predictions
cts_bj <-
  "out/BaseJumper/bj-expression/PD63118_250224_152704/tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv" %>%
  readr::read_tsv() %>%
  dplyr::select(id = SampleId, celltype_bj = Progenitor) %>%
  tidyr::pivot_longer(cols = "celltype_bj")

# load immunelens celltype predictions
scores <-
  readr::read_tsv("out/immunelens/scores.tsv") %>%
  dplyr::filter(gene != "TCRB")

# plot scores
p_il <-
  scores %>%
  ggplot(aes(x = id, y = gene, fill = cell_fraction)) +
  geom_tile() +
  guides(x = guide_axis(angle = -90))

# combine predictions
cts <-
  dplyr::bind_rows(cts_vdj, cts_bj) %>%
  dplyr::full_join(ss_exp) %>%
  dplyr::filter(name %in% c("celltype_VDJ_recomb", "celltype_SHM") & seq_type == "dna" |
                name == "celltype_bj" & seq_type == "rna") %>%
  # harmonise celltype names
  dplyr::mutate(
    ct = dplyr::case_when(
      value %in% c("B cell", "B_cell", "B mem") ~ "B cell",
      value %in% c("T cell", "T_cells", "alpha-beta T cell") ~ "T cell",
      value %in% c("Chondrocytes", "Fibroblasts", "Smooth_muscle_cells") ~ "other",
      TRUE ~ value
    ))

# plot
p <-
  cts %>%
  ggplot(aes(x = cell_id, y = name, fill = ct)) +
  geom_tile() +
  # rotate x axis labels
  guides(x =  guide_axis(angle = -90))

# plot ct props
cts %>%
  dplyr::mutate(ct = ifelse(is.na(ct), "unannotated", ct)) %>%
  ggplot(aes(x = ct, fill = ct)) +
  geom_bar() +
  ggh4x::facet_grid2(~ name, scales = "free_x", space = "free_x")
