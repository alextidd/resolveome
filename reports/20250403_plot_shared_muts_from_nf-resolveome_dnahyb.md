---
title: "Shared mutations from the nf-resolveome genotyping output"
author: "Alexandra Tidd"
date: "03 April, 2025"
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
data_dir <- "out/nf-resolveome/dnahyb/"
man_insp <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  filter((!chr_dropout | is.na(chr_dropout)) & (!suspected_doublet | is.na(suspected_doublet))) %>%
  mutate(loh_1p = as.character(as.numeric(`1p_loh`))) %>%
  distinct(cell_id, loh_1p, celltype_SHM, celltype_VDJ_recomb, TNFRSF14_mut)
ss <-
  readr::read_csv(file.path(data_dir, "samplesheet.csv")) %>%
  left_join(man_insp)
```

We load the variants.


``` r
genos <-
  readr::read_tsv(file.path(data_dir, "PD63118/genotyping/mutations/PD63118_annotated_mutations.tsv")) %>%
  left_join(ss)
```


We generate a heatmap of all mutations.


``` r
p <- plot_vaf_heatmap(genos, show_rownames = FALSE, p_title = "all mutations")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250403_plot_shared_muts_from_nf-resolveome_dnahyb_files/figure-html/plot_heatmap-1.png)<!-- -->

We generate a heatmap of the shared mutations.


``` r
p <-
  genos %>%
  filter(mut_depth > 0) %>%
  group_by(chr, pos, ref, alt) %>%
  filter(n() > 1) %>%
  plot_vaf_heatmap(show_rownames = FALSE, p_title = "shared mutations")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250403_plot_shared_muts_from_nf-resolveome_dnahyb_files/figure-html/plot_heatmap_shared-1.png)<!-- -->

We filter to `mut_vaf` > 0.25 and `mut_depth` > 1.


``` r
genos_filter <-
  genos %>%
  filter(mut_vaf > 0.2, mut_depth > 1)

# return
genos_filter
```

```
## # A tibble: 158 × 40
##    donor_id chr       pos ref       alt   gene  strand ref_cod mut_cod ref3_cod mut3_cod aachange ntchange codonsub impact pid   total_depth ref_depth
##    <chr>    <chr>   <dbl> <chr>     <chr> <chr>  <dbl> <chr>   <chr>   <chr>    <chr>    <chr>    <chr>    <chr>    <chr>  <chr>       <dbl>     <dbl>
##  1 PD63118  1     2488098 TGAGGCAT… T     TNFR…      1 .       .       .        .        .        1-8-del… .        no-SNV ENSP…          18        11
##  2 PD63118  1     2488105 T         A     TNFR…      1 T       A       ATG      AAG      M1K      T2A      ATG>AAG  Misse… ENSP…         131         0
##  3 PD63118  1     2488106 G         T     TNFR…      1 G       T       TGG      TTG      M1I      G3T      ATG>ATT  Misse… ENSP…         137        78
##  4 PD63118  1     2488138 G         A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP…         104        59
##  5 PD63118  1     2488138 G         A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP…          28         3
##  6 PD63118  1     2488138 G         A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP…         124        48
##  7 PD63118  1     2488138 G         A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP…         157        98
##  8 PD63118  1     2488139 G         A     TNFR…      1 G       A       GGA      GAA      W12*     G36A     TGG>TGA  Nonse… ENSP…          87        51
##  9 PD63118  1     2488139 G         A     TNFR…      1 G       A       GGA      GAA      W12*     G36A     TGG>TGA  Nonse… ENSP…         122        70
## 10 PD63118  1     2488139 G         A     TNFR…      1 G       A       GGA      GAA      W12*     G36A     TGG>TGA  Nonse… ENSP…          94        24
## # ℹ 148 more rows
## # ℹ 22 more variables: mut_depth <dbl>, mut_vaf <dbl>, type <chr>, id <chr>, cell_id <chr>, plate <dbl>, well <chr>, seq_type <chr>, run_id <dbl>,
## #   lane <chr>, plex_n <dbl>, study_id <dbl>, sanger_sample_id <chr>, supplier_sample_name <chr>, manifest_file <chr>, bam <chr>, mutations <chr>,
## #   snps <chr>, loh_1p <chr>, celltype_SHM <chr>, celltype_VDJ_recomb <chr>, TNFRSF14_mut <chr>
```


``` r
p <- plot_vaf_heatmap(genos_filter, p_title = "all mutations")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250403_plot_shared_muts_from_nf-resolveome_dnahyb_files/figure-html/plot_heatmap_filtered-1.png)<!-- -->


``` r
p <- genos_filter %>%
  group_by(chr, pos, ref, alt) %>%
  filter(n() > 1) %>%
  plot_vaf_heatmap(p_title = "shared mutations")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250403_plot_shared_muts_from_nf-resolveome_dnahyb_files/figure-html/plot_heatmap_filtered_shared-1.png)<!-- -->

We generate a heatmap of the nonsynonymous mutations.


``` r
p <-
  genos_filter %>%
  filter(impact != "Synonymous") %>%
  plot_vaf_heatmap("non-synonymous mutations", annotations = "loh_1p")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250403_plot_shared_muts_from_nf-resolveome_dnahyb_files/figure-html/plot_nonsym-1.png)<!-- -->

We generate a heatmap of the TNFRSF14 mutations.


``` r
p <-
  genos_filter %>%
  filter(gene == "TNFRSF14") %>%
  plot_vaf_heatmap("TNFRSF14 mutations")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250403_plot_shared_muts_from_nf-resolveome_dnahyb_files/figure-html/plot_heatmap_filtered_TNFRSF14-1.png)<!-- -->
