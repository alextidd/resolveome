---
title: "ResolveOME celltype prediction"
author: "Alexandra Tidd"
date: "02 April, 2025"
output:
  html_document:
    fig_width: 8
    keep_md: true
    toc: true
    toc_float: true
    toc_collapsed: true
    toc_depth: 4
    theme: lumen
params:
  out_dir: null
---



First, we load the samplesheet with all metadata for the cells.


``` r
# load samplesheet
ss <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::mutate(run_id = as.character(run_id))
ss_exp <-
  ss %>%
  dplyr::select(id, cell_id, seq_type) %>%
  dplyr::group_by(id, cell_id) %>%
  tidyr::expand_grid(name = c("celltype_SHM", "celltype_VDJ_recomb", "celltype_bj_expression")) %>%
  dplyr::ungroup()
```

Next, we load the celltype predictions from BaseJumper, which uses SingleR to
make predictions based on the Human Primary Cell Atlas (HPCA).


``` r
# get latest bj dir
bj_dir <-
  list.files(params$out_dir, pattern = "PD63118_",
             include.dirs = TRUE, full.names = TRUE) %>%
  sort() %>%
  tail(1)

# get basejumper celltype preds
cts_bj <-
  readr::read_tsv(paste0(bj_dir, "/tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv")) %>%
  dplyr::left_join(ss, by = c("SampleId" = "id"))
```

## BaseJumper bj-expression predictions


``` r
# plot celltype counts
cts_bj %>%
  dplyr::mutate(plate = as.character(plate)) %>%
  dplyr::add_count(Progenitor, name = "total") %>%
  dplyr::count(Progenitor, plate, total) %>%
  ggplot(aes(x = reorder(Progenitor, -total), y = n, fill = plate)) +
  geom_col() +
  geom_text(aes(y = total, label = total), vjust = -0.5) +
  guides(x = guide_axis(angle = -90)) +
  theme_classic() +
  labs(x = "SingleR HPCA celltype prediction") +
  scale_fill_brewer(palette = "Set1")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250402_bj-expression_celltype_prediction_files/figure-html/plot_bj-1.png)<!-- -->
