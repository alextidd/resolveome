---
title: "Shared mutations from the nf-resolveome genotyping output"
author: "Alexandra Tidd"
date: "29 March, 2025"
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



## DNA

We load the samplesheet.


``` r
data_dir <- "out/nf-resolveome/dna/"
man_insp <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::filter(plate == 3) %>%
  dplyr::distinct(cell_id, `1p_loh`, celltype_SHM, celltype_VDJ_recomb)
ss <-
  readr::read_csv(file.path(data_dir, "samplesheet.csv")) %>%
  dplyr::left_join(man_insp)
```

We load the variants.


``` r
geno <-
  readr::read_tsv(file.path(data_dir, "PD63118/genotyping/mutations/PD63118_annotated_mutations.tsv")) %>%
  dplyr::left_join(ss)
```

We generate a heatmap of the shared mutations.


``` r
# plot shared mutations
geno %>%
  plot_vaf_heatmap("NanoSeq mutations - DNA - all")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250319_plot_shared_muts_from_nf-resolveome_dna_files/figure-html/plot_heatmap-1.png)<!-- -->



``` r
# plot shared mutations
geno %>%
  dplyr::filter(mut_depth > 0) %>%
  dplyr::group_by(chr, pos, ref, alt) %>%
  dplyr::filter(dplyr::n_distinct(id) > 1) %>%
  dplyr::ungroup() %>%
  plot_vaf_heatmap("NanoSeq mutations - DNA - shared")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250319_plot_shared_muts_from_nf-resolveome_dna_files/figure-html/plot_heatmap_shared-1.png)<!-- -->

## DNA hyb

We load the samplesheet.


``` r
ss <-
  readr::read_csv("out/nf-resolveome/dnahyb/samplesheet.csv") %>%
  dplyr::left_join(man_insp)
```

We load the variants.


``` r
geno_hyb <-
  readr::read_tsv("out/nf-resolveome/dnahyb/PD63118/genotyping/mutations/PD63118_annotated_mutations.tsv") %>%
  dplyr::left_join(ss)
```

We plot a histogram of the VAFs of the mutations.


``` r
geno_hyb %>%
  dplyr::filter(mut_depth > 0) %>%
  ggplot(aes(x = mut_vaf)) +
  geom_histogram(binwidth = 0.05) +
  theme_classic()
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250319_plot_shared_muts_from_nf-resolveome_dna_files/figure-html/plot_vaf-1.png)<!-- -->

``` r
geno_hyb %>%
  dplyr::arrange(total_depth) %>%
  dplyr::count(total_depth) %>%
  dplyr::mutate(prop_cov = n / dplyr::n()) %>%
  ggplot(aes(x = total_depth, y = prop_cov)) +
  geom_step()
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250319_plot_shared_muts_from_nf-resolveome_dna_files/figure-html/plot_vaf-2.png)<!-- -->

We generate a heatmap of the mutations.


``` r
# plot shared mutations
geno_hyb %>%
  dplyr::filter(mut_depth > 0) %>%
  plot_vaf_heatmap("NanoSeq mutations - DNAhyb - all")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250319_plot_shared_muts_from_nf-resolveome_dna_files/figure-html/plot_heatmap_hyb-1.png)<!-- -->


``` r
# plot shared mutations
geno_hyb %>%
  dplyr::filter(mut_depth > 0) %>%
  dplyr::group_by(chr, pos, ref, alt) %>%
  dplyr::filter(dplyr::n_distinct(id) > 1) %>%
  dplyr::ungroup() %>%
  plot_vaf_heatmap("NanoSeq mutations - DNAhyb - shared")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250319_plot_shared_muts_from_nf-resolveome_dna_files/figure-html/plot_heatmap_hyb_shared-1.png)<!-- -->
