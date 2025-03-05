---
title: "Shared mutations from the bj-somatic-variantcalling WGS output"
author: "Alexandra Tidd"
date: "05 March, 2025"
output:
  html_document:
    fig_width: 8
    keep_md: true
    toc: true
    toc_float: true
    toc_collapsed: true
    toc_depth: 4
    theme: lumen
---



We load the samplesheet.


``` r
man_insp <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::mutate(run_id = as.character(run_id))
metadata <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::left_join(man_insp)
ss <-
  readr::read_csv("out/BaseJumper/bj-somatic-variantcalling/wgs/samplesheet.csv")
```

We load the annotated variants.


``` r
# get the latest run
bj_dir <-
  list.files("out/BaseJumper/bj-somatic-variantcalling/wgs/",
             pattern = "PD63118_", include.dirs = TRUE, full.names = TRUE) %>%
  sort() %>%
  tail(2) %>% head(1)
vcfs_subdir <- "/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SUBSET_VCF_VARIANTS/"

# get vafs
vafs <-
  ss$biosampleName %>%
  purrr::set_names() %>%
  purrr::map(function(id) {
    vcf_file <- paste0(bj_dir, vcfs_subdir, id, "_somatic_filtered_variants.vcf.gz")
    if (file.exists(vcf_file)) {
      vcf <-
        readr::read_tsv(vcf_file, comment = "##") %>%
        {split(., .$FORMAT)} %>%
        purrr::map(function(vcf_i) {
          format_i <- stringr::str_split_1(unique(vcf_i$FORMAT), ":")
          vcf_i %>%
            dplyr::mutate(gt = get(id)) %>%
            tidyr::separate_wider_delim("gt", delim = ":", names = format_i) %>%
            janitor::clean_names() %>%
            dplyr::mutate(chr = number_chrom, total_depth = dp,
                          allele = paste(ref, alt, sep = ","),
                          alt = gsub(",.*", "", alt)) %>%
            tidyr::separate_longer_delim(cols = c("allele", "ad"), delim = ",") %>%
            dplyr::select(chr, pos, ref, alt, allele, total_depth, ad) %>%
            dplyr::mutate(name = dplyr::case_when(allele == ref ~ "ref_depth",
                                                  allele == alt ~ "alt_depth",
                                                  TRUE ~ "other_depth")) %>%
            dplyr::filter(name != "other_depth") %>%
            tidyr::pivot_wider(
              names_from = "name", values_from = "ad",
              id_cols = c("chr", "pos", "ref", "alt", "total_depth")) %>%
            readr::type_convert() %>%
            dplyr::mutate(alt_vaf = alt_depth / total_depth)
        }) %>%
        dplyr::bind_rows()
    }
  }) %>%
  purrr::compact() %>%
  dplyr::bind_rows(.id = "id")

# save vafs
vafs %>% readr::write_tsv("out/analysis/wgs_bj-somatic-variantcalling_mut_vafs.tsv")
```

We remove cells with chromosomal dropout and suspected doublets.


``` r
vafs <-
  vafs %>%
  dplyr::left_join(metadata) %>%
  dplyr::filter(!chr_dropout, !suspected_doublet)
```

We annotate them with `dndscv`.


``` r
vafs <-
  vafs %>%
  dplyr::transmute(sampleID = id, chr = gsub("chr", "", chr),
                   pos, ref, mut = alt) %>%
  dplyr::distinct() %>%
  dndscv(refdb = "../reference/dndscv/RefCDS_human_GRCh38_GencodeV18_recommended.rda",
         max_coding_muts_per_sample = Inf,
         max_muts_per_gene_per_sample = Inf,
         outp = 1) %>%
  {.$annotmuts} %>%
  dplyr::transmute(chr = paste0("chr", chr), pos, ref, alt = mut, gene,
                   id = sampleID, aachange, ntchange, codonsub, pid, impact) %>%
  dplyr::full_join(vafs) %>%
  dplyr::mutate(impact = ifelse(is.na(impact), "Non-coding", impact) %>%
                  factor(levels = rev(names(impact_colours))))
```

Plot the distribution of mutation types per cell.


``` r
# all muts
vafs %>%
  dplyr::add_count(id) %>%
  ggplot(aes(x = reorder(id, -n), fill = impact)) +
  geom_bar() +
  scale_fill_manual(values = impact_colours) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = guide_axis(angle = -90)) +
  labs(x = "")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_mut_types-1.png)<!-- -->

``` r
# coding muts
vafs %>%
  dplyr::filter(impact != "Non-coding") %>%
  dplyr::add_count(id) %>%
  ggplot(aes(x = reorder(id, -n), fill = impact)) +
  geom_bar() +
  scale_fill_manual(values = impact_colours) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = guide_axis(angle = -90)) +
  labs(x = "")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_mut_types-2.png)<!-- -->

Plot the VAF distribution.


``` r
c(1, 5, 10, 20, 50) %>%
  purrr::walk(function(min_total_depth) {
    p_dat <- vafs %>% dplyr::filter(total_depth >= min_total_depth)
    p <-
      p_dat %>%
      ggplot(aes(x = alt_vaf)) +
      geom_histogram() +
      theme_classic() +
      labs(title = paste("VAF distribution, min depth =", min_total_depth),
           subtitle = paste(nrow(p_dat), "mutations"), x = "alt VAF")
    print(p)
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf-3.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf-4.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf-5.png)<!-- -->

Plot a heatmap of the shared variants.


``` r
# function: plot VAF heatmap
plot_vaf_heatmap <- function(p_dat, p_title, show_rownames = FALSE) {
  # reshape data for heatmap
  heatmap_data <-
    p_dat %>%
    reshape2::dcast(mut_id + n_cells_w_mut ~ id, value.var = "alt_vaf")
  rownames(heatmap_data) <- heatmap_data$mut_id
  heatmap_matrix <- as.matrix(heatmap_data[, -c(1:3)]) # drop annotations
  n_total_muts <- dplyr::n_distinct(p_dat$mut_id)
  n_total_cells <- dplyr::n_distinct(p_dat$id)

  # calculate the number of mutations per column (id)
  annotation_col <-
    p_dat %>%
    dplyr::count(id, name = "n_muts") %>%
    dplyr::left_join(metadata %>% dplyr::select(run_id, celltype_VDJ_recomb, id)) %>%
    dplyr::mutate(run_id = as.character(run_id)) %>%
    tibble::column_to_rownames("id")

  # replace NAs with 0
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # prepare annotations for columns
  annotation_row <- data.frame(n = heatmap_data$n)
  rownames(annotation_row) <- heatmap_data$mut_id

  # plot
  pheatmap::pheatmap(
    mat = heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    annotation_row = annotation_row,
    annotation_col = annotation_col,
    color = colorRampPalette(c("white", "blue"))(50),
    main = paste("VAF heatmap:", p_title, "\n",
                 n_total_cells, "cells,", n_total_muts, "mutations"))
}
```


``` r
# plot all shared muts
p_dat_shared <-
  vafs %>%
  dplyr::mutate(mut_id = paste(chr, pos, ref, alt, sep = "-"),
                n_cells = dplyr::n_distinct(id)) %>%
  dplyr::add_count(mut_id, name = "n_cells_w_mut") %>%
  dplyr::mutate(prop_cells_w_mut = n_cells_w_mut / n_cells) %>%
  # shared
  dplyr::filter(n_cells_w_mut > 1)
p <- plot_vaf_heatmap(p_dat_shared, "all shared muts")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_shared_vaf_heatmap-1.png)<!-- -->


``` r
# look at shared muts, exclude those widely shared
p_dat_shared_filtered <-
  p_dat_shared %>%
  dplyr::filter(prop_cells_w_mut < 0.3)
p <- plot_vaf_heatmap(p_dat_shared_filtered, "shared muts, % cells w mut < 0.3")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_shared_vaf_heatmap_filtered-1.png)<!-- -->


``` r
# look at coding muts
p_dat_coding <-
  p_dat_shared %>%
  dplyr::filter(impact != "Non-coding")
p <- plot_vaf_heatmap(p_dat_coding, "shared coding muts", TRUE)
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_shared_vaf_heatmap_coding-1.png)<!-- -->
