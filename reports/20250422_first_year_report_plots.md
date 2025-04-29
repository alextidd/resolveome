---
title: "First year report plots"
author: "Alexandra Tidd"
date: "28 April, 2025"
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
  resume: true
---





# Samplesheet


``` r
man_insp <-
  read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::select(-run_id) %>%
  # fix loh_1p
  mutate(loh_1p = as.character(as.numeric(loh_1p)))
ss <-
  read_csv("data/resolveome/samplesheet_local.csv") %>%
  # add manual inspection data
  left_join(man_insp)
driver_genes <-
  readLines("../trencadis-seq/data/thyroid/driver_genes/driver_genes.txt")
```

# Quality control

We determine which cells to keep. QC has been completed for plate 3, but has not
been completed for plate 10, as only DNAhyb data is available for these cells,
so it is more difficult to definitively call doublets and chromosomal dropouts.
However, based on the new gating strategy for flow cytometry and based on the
pre-PCR quants, we can be fairly confident that the cells are of high quality,
so we will include all plate 10 cells in the analysis.


``` r
# counts
n_seq_cells <- n_distinct(ss$cell_id)
n_seq_dnahyb_cells <- ss %>% filter(seq_type == "dnahyb") %>% {n_distinct(.$cell_id)}
n_seq_dna_cells <- ss %>% filter(seq_type == "dna") %>% {n_distinct(.$cell_id)}

# qc-passing cells
ss_pass <-
  ss %>%
  filter((plate == 3 & !suspected_doublet & !chr_dropout) |
         (plate == 10 &
          (is.na(suspected_doublet) | !suspected_doublet) &
          (is.na(chr_dropout) | !chr_dropout)))
n_pass_cells <- n_distinct(ss_pass$cell_id)
```

Of 263 isolated nuclei, 182 passed pre-PCR quantification QC
and were submitted for sequencing. 182 cells underwent 
targeted DNA and RNA sequencing. 99/182 cells
additionally underwent whole genome DNA sequencing. 115 cells 
passed QC.

# Celltype annotation

## VDJ + SHM + CSR


``` r
cts <-
  ss_pass %>%
  # synthesise celltype annotations from VDJ + SHM + CSR
  mutate(
    celltype_VDJ_SHM_CSR = case_when(
      celltype_VDJ_recomb == "B cell" &
        (celltype_SHM == "mature B cell" |
         class_switch_recombination_CSR == TRUE) ~ "activated B",
      celltype_VDJ_recomb == "B cell" & celltype_SHM == "not mature B cell" ~ "naive B",
      celltype_VDJ_recomb == "alpha-beta T cell" ~ "T",
      celltype_VDJ_recomb == "not lymphocyte" ~ "other",
      TRUE ~ NA) %>%
      factor(levels = c("activated B", "naive B", "T", "other"))) %>%
  filter(!is.na(celltype_VDJ_SHM_CSR)) %>%
  distinct(cell_id, celltype_VDJ_SHM_CSR)

# subset samplesheet
ss_annot <-
  ss_pass %>%
  inner_join(cts)

# plot
cts %>%
  count(celltype_VDJ_SHM_CSR) %>%
  ggplot(aes(x = reorder(celltype_VDJ_SHM_CSR, -n), y = n, fill = celltype_VDJ_SHM_CSR)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5) +
  scale_fill_manual(values = pal) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = "VDJ + SHM celltype")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/get_celltypes-1.png)<!-- -->

## ImmuneLENS


``` r
# immune lens annotations
immunelens <-
  "out/immunelens/scores.tsv" %>%
  read_tsv() %>%
  inner_join(ss_annot) %>%
  transmute(cell_id, name = paste0("cell_frac_", gene), cell_fraction,
            celltype_VDJ_SHM_CSR)

# plot
immunelens %>%
  filter(name %in% c("cell_frac_TCRA", "cell_frac_IGH")) %>%
  ggplot(aes(x = celltype_VDJ_SHM_CSR, y = cell_fraction, colour = celltype_VDJ_SHM_CSR)) +
  geom_boxplot() +
  geom_jitter(height = 0) +
  facet_grid(~ name) +
  scale_colour_manual(values = pal) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/il_annot-1.png)<!-- -->

## bj-expression


``` r
bje <-
  "out/BaseJumper/bj-expression/rna/PD63118_run/tertiary_analyses/classification_cell_typing/df_cell_typing_summary_singler_hpca_gtex_tcga.tsv" %>%
  read_tsv() %>%
  dplyr::rename(id = SampleId) %>%
  inner_join(ss_annot)

bje %>% count(Progenitor)
```

```
## # A tibble: 15 × 2
##    Progenitor               n
##    <chr>                <int>
##  1 Astrocyte                3
##  2 B_cell                  48
##  3 Chondrocytes             2
##  4 Epithelial_cells         1
##  5 Fibroblasts              1
##  6 Hepatocytes              1
##  7 Macrophage               1
##  8 Monocyte                 1
##  9 NK_cell                  4
## 10 Neuroepithelial_cell     1
## 11 Neurons                  3
## 12 Osteoblasts              6
## 13 Platelets                6
## 14 Smooth_muscle_cells      1
## 15 T_cells                 27
```

``` r
bje %>% filter(celltype_VDJ_SHM_CSR == "other") %>% count(Progenitor)
```

```
## # A tibble: 6 × 2
##   Progenitor           n
##   <chr>            <int>
## 1 Epithelial_cells     1
## 2 Macrophage           1
## 3 Monocyte             1
## 4 NK_cell              1
## 5 Osteoblasts          2
## 6 Platelets            2
```

``` r
bje %>%
  mutate(celltype_bje = case_when(
    Progenitor %in% c("B_cell", "T_cells", "NK_cell", "Macrophage", "Monocyte") ~ Progenitor,
    TRUE ~ "other")) %>%
  count(celltype_VDJ_SHM_CSR, celltype_bje)
```

```
## # A tibble: 14 × 3
##    celltype_VDJ_SHM_CSR celltype_bje     n
##    <fct>                <chr>        <int>
##  1 activated B          B_cell          43
##  2 activated B          NK_cell          1
##  3 activated B          T_cells          1
##  4 activated B          other           11
##  5 naive B              B_cell           2
##  6 naive B              other            2
##  7 T                    B_cell           3
##  8 T                    NK_cell          2
##  9 T                    T_cells         26
## 10 T                    other            7
## 11 other                Macrophage       1
## 12 other                Monocyte         1
## 13 other                NK_cell          1
## 14 other                other            5
```


``` r
ct_stats <-
  cts %>% count(celltype_VDJ_SHM_CSR) %>% {split(.$n, .$celltype_VDJ_SHM_CSR)}

il_stats <-
  immunelens %>%
  mutate(
    celltype = case_when(celltype_VDJ_SHM_CSR %in% c("activated B", "naive B") ~ "B",
                         TRUE ~ celltype_VDJ_SHM_CSR)) %>%
  group_by(celltype, name) %>%
  summarise(mean = mean(cell_fraction, na.rm = TRUE),
            max = max(cell_fraction, na.rm = TRUE),
            min = min(cell_fraction, na.rm = TRUE)) %>%
  mutate(stat = paste0(round(mean * 100), "% (", round(min * 100), "-",
                       round(max * 100), "%)")) %>%
  split(.$name) %>%
  purrr::map(~ split(.x$stat, .x$celltype))

bje_stats <-
  bje %>%
  left_join(cts) %>%
  transmute(cell_id, Progenitor,
            celltype_bje = case_when(Progenitor == "B_cell" ~ "B",
                                     Progenitor == "T_cells" ~ "T",
                                     TRUE ~ "other"),
            celltype = case_when(celltype_VDJ_SHM_CSR %in% c("activated B", "naive B") ~ "B",
                                 TRUE ~ celltype_VDJ_SHM_CSR),
            lineage = ifelse(celltype == "other", "other", "B/T")) %>%
  count(lineage, celltype, celltype_bje)
bje_stats %>% knitr::kable()
```



|lineage |celltype |celltype_bje |  n|
|:-------|:--------|:------------|--:|
|B/T     |B        |B            | 45|
|B/T     |B        |T            |  1|
|B/T     |B        |other        | 14|
|B/T     |T        |B            |  3|
|B/T     |T        |T            | 26|
|B/T     |T        |other        |  9|
|other   |other    |other        |  8|

``` r
bje_non_BT <-
  bje %>%
  left_join(cts) %>%
  transmute(cell_id, Progenitor, celltype_VDJ_SHM_CSR) %>%
  filter(celltype_VDJ_SHM_CSR == "other") %>%
  group_by(Progenitor) %>%
  summarise(n = n(), cell_ids = paste(cell_id, collapse = ",")) %>%
  arrange(-n)
bje_non_BT %>% knitr::kable()
```



|Progenitor       |  n|cell_ids                      |
|:----------------|--:|:-----------------------------|
|Osteoblasts      |  2|plate3_wellA9,plate3_wellD12  |
|Platelets        |  2|plate3_wellB10,plate3_wellB12 |
|Epithelial_cells |  1|plate10_wellC10               |
|Macrophage       |  1|plate3_wellA4                 |
|Monocyte         |  1|plate3_wellA5                 |
|NK_cell          |  1|plate10_wellD5                |

``` r
n_bt <- bje_stats %>% filter(lineage == "B/T") %>% pull(n) %>% sum()
n_bt_correct <- bje_stats %>% filter(lineage == "B/T" & celltype == celltype_bje) %>% pull(n) %>% sum()
n_other <- bje_stats %>% filter(lineage == "other") %>% pull(n) %>% sum()
```

Based on VDJ, SHM, and CSR annotations, 
60 cells were B cells 
(56 of which were activated), 
38 were T cells, and 8 were other celltypes. 

`ImmuneLENS` scores for T cell fraction based on the TCRA locus strongly confirm 
the T cell annotations with a mean of 96% (95-98%), versus a
mean of 10% (0-41%) for B cells.
Mean ImmuneLENS scores based on the IGH were 79% (54-96%) for B
cells and 38% (22-59%) for T cells. This correlated well with B
cell annotation, although did not differentiate B and T cells as effectively, 
with some overlap in ranges [see supplementary figure X]. 

`bj-expression` automated labelling of the RNA data was compared to the
DNA-based annotations and was found to perform quite sporadically. 
71/98 B and T lymphocytes were annotated correctly.
Among the 8 cells that were not B or T cells according to the DNA
data, there were 2 platelets, 1 macrophage, 1 monocyte, 1 NK cell, and 1 
epithelial cell, which are plausible labels. However, two non-B/T cells were 
labelled as osteoblasts, which is an unrealistic label. Nonetheless, most of
these automated labels show consensus with the DNA-based annotations.

Final celltypes were assigned using the outcomes of the VDJ, SHM, and CSR
inspections, which were largely supported by the `ImmuneLENS` and
`bj-expression` results.

# Sequencing metrics

We load the DNA sequencing coverage.


``` r
# get paths to mosdepth output
mosdepth <-
  ss_annot %>%
  transmute(
    seq_type, donor_id, id,
    out_dir = file.path("out/nf-resolveome", seq_type, donor_id, id),
    summary_txt = paste0(out_dir, "/mosdepth/", id, ".mosdepth.summary.txt"),
    dist = paste0(out_dir, "/mosdepth/", id, ".mosdepth.",
                  ifelse(seq_type == "dna", "global", "region"),
                  ".dist.txt")) %>%
  filter(file.exists(summary_txt), file.exists(dist))

# load mean coverage
cov <-
  xfun::cache_rds({
    cov <-
      split(mosdepth$summary_txt, mosdepth$id) %>%
      purrr::map(read_tsv) %>%
      bind_rows(.id = "id") %>%
      dplyr::rename(chr = chrom) %>%
      filter(!grepl("^GL", chr), !grepl("Un", chr), !grepl("random|HLA", chr),
            !grepl("alt$", chr), chr %in% c("total", "total_region")) %>%
      left_join(ss %>% dplyr::select(id, seq_type)) %>%
      filter(chr == "total" & seq_type == "dna" |
            chr == "total_region" & seq_type == "dnahyb") %>%
      transmute(id, name = "mean coverage", value = mean)
    cov
  }, file = "cov.rds", resume = params$resume)

# load coverage distribution
dist <-
  split(mosdepth$dist, mosdepth$id) %>%
  purrr::map(read_tsv, col_names = c("chr", "cov", "prop")) %>%
  bind_rows(.id = "id") %>%
  filter(chr == "total") %>%
  inner_join(ss_annot)
prop_cov <-
  dist %>%
  filter(cov == 1) %>%
  transmute(id, name = "prop of genome covered", value = prop)
```

We load the RNA sequencing counts and features.


``` r
# dirs
bje_dir <- "out/BaseJumper/bj-expression/rna/PD63118_run/"

# load gene counts
gene_counts <-
  file.path(bje_dir, "secondary_analyses/quantification_salmon/df_gene_counts_salmon.tsv") %>%
  read_tsv() %>%
  mutate(id = gsub("salmon_outdir_", "", File)) %>%
  group_by(id) %>%
  summarise(value = n_distinct(gene_symbol_gene_id), name = "n genes")

# get read counts
read_counts <-
  file.path(bje_dir, "secondary_analyses/quantification_salmon/df_gene_counts_salmon.tsv") %>%
  read_tsv() %>%
  mutate(id = gsub("salmon_outdir_", "", File)) %>%
  group_by(id) %>%
  summarise(value = sum(countsFromAbundanceNo), name = "n counts")
```

We plot the coverage distribution.


``` r
dist %>%
  filter(prop >= 0.01) %>%
  ggplot(aes(x = cov, y = prop, colour = cell_id)) +
  geom_line() +
  theme_classic() +
  labs(x = "coverage", y = "proportion of genome") +
  facet_grid(~ seq_type, scales = "free_x") +
  theme(legend.position = "none")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_coverage_dist-1.png)<!-- -->

We plot the coverage.


``` r
p_dat <-
  bind_rows(read_counts, gene_counts, cov, prop_cov) %>%
  inner_join(ss_annot) %>%
  group_by(cell_id) %>%
  mutate(seq = paste(sort(unique(seq_type)), collapse = " + ")) %>%
  group_by(seq) %>%
  mutate(seq_n = paste0(seq, "\n(n = ", n_distinct(cell_id), ")"))

p_dat %>%
  left_join(
    bind_rows(
      p_dat %>%
        filter(name == "mean coverage", seq == "dna + dnahyb + rna") %>%
        distinct(cell_id, cov = value),
      p_dat %>%
        filter(name == "mean coverage", seq == "dnahyb + rna") %>%
        distinct(cell_id, cov = value)
    )) %>%
  ggplot(aes(x = reorder(cell_id, -cov), y = value)) +
  geom_col(position = "dodge") +
  ggh4x::facet_nested(seq_type + name ~ seq_n, scales = "free", space = "free_x", switch = "y") +
  guides(x = guide_axis(angle = -90)) +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.placement = "outside",
        axis.title.y = element_blank(), axis.title.x = element_blank())
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_coverage-1.png)<!-- -->


``` r
seq_stats <-
  p_dat %>%
  group_by(name, seq_type) %>%
  summarise(mean = mean(value), min = min(value), max = max(value)) %>%
  pivot_longer(cols = c(mean, min, max), names_to = "stat",
               values_to = "value") %>%
  split(.$seq_type) %>%
  purrr::map(~ split(.x, .x$name) %>% purrr::map(~ split(.x$value, .x$stat)))
```

Mean proportion coverage was
92%
for DNA genome-wide, and
95%
for targeted DNA within the immune panel.
Average coverage was
25X (range
18X -
35X) 
for the DNA genome-wide, and
324X (range
37X -
518X)
for the targeted DNA within the immune panel [see supplementary figure X].

In the RNA sequencing, the average number of unique genes counted per cell was
2,658 (range
1,123 -
12,001)
and the average total count per cell was
2,342,321 (range
32,661 -
16,044,163)
[see supplementary figure X].

# Copy number and LOH profiles


``` r
loh_counts <-
  ss_annot %>%
  distinct(cell_id, loh_1p) %>%
  count(loh_1p) %>%
  {split(.$n, .$loh_1p)}

dna_loh_counts <-
  ss_annot %>%
  filter(seq_type == "dna") %>%
  distinct(cell_id, loh_1p) %>%
  count(loh_1p) %>%
  {split(.$n, .$loh_1p)}

ss_annot %>%
  distinct(cell_id, loh_1p, celltype_VDJ_SHM_CSR) %>%
  count(loh_1p, celltype_VDJ_SHM_CSR) %>%
  knitr::kable()
```



|loh_1p |celltype_VDJ_SHM_CSR |  n|
|:------|:--------------------|--:|
|0      |activated B          | 29|
|0      |naive B              |  4|
|0      |T                    | 38|
|0      |other                |  8|
|1      |activated B          | 27|

27 / 106 cells had
evidence of 1p LOH. All cells with 1p LOH were activated B cells.

For cells with whole genome DNA sequencing, there is sufficient resolution in 
the BAF plots to manually place the breakpoints. Of the 
39
cells with whole genome DNA sequencing,
11 
cells had evidence of 1p LOH. At least 8 distinct 1p
breakpoints were confirmed among these 
11 
cells, suggesting 11+ independent 1p 
copy number events in a population of just 
39
cells [see figure X]. The Ginkgo copy number calls show no signal of 1p coverage 
loss across these cells [see supplementary figure X], suggesting that these are 
all copy neutral LOH events. 

# Somatic mutation genotyping


``` r
# load the genotyping and pile up reads across DNA and targeted DNA.
genos <-
  ss %>%
  # get genotyping for dna and dnahyb
  filter(seq_type %in% c("dna", "dnahyb")) %>%
  distinct(seq_type, donor_id) %>%
  mutate(genos = paste0("out/nf-resolveome/", seq_type, "/", donor_id,
                        "/genotyping/mutations/", donor_id,
                        "_annotated_mutations.tsv")) %>%
  purrr::pmap(function(seq_type, donor_id, genos) {
    read_tsv(genos) %>%
      mutate(donor_id = donor_id) %>%
      # convert id -> cell_id
      inner_join(ss %>% dplyr::select(id, cell_id)) %>%
      dplyr::select(-id)
  }) %>%
  # pile up
  bind_rows() %>%
  dplyr::select(-mut_vaf) %>%
  group_by(across(-matches("_depth$"))) %>%
  summarise(across(matches("_depth$"), \(x) sum(x, na.rm = TRUE))) %>%
  ungroup() %>%
  mutate(mut_vaf = mut_depth / total_depth) %>%
  # annotate mutations on 1p
  mutate(on_1p = ifelse(chr == "1" & pos < 123400000, TRUE, FALSE)) %>%
  # filter to only annotated, qc'ed cells
  inner_join(ss_annot %>% distinct(cell_id, celltype_VDJ_SHM_CSR, loh_1p))

# return
genos
```

```
## # A tibble: 152,004 × 25
##    donor_id chr       pos ref   alt   gene  strand ref_cod mut_cod ref3_cod mut3_cod aachange ntchange codonsub impact pid   type  cell_id total_depth
##    <chr>    <chr>   <dbl> <chr> <chr> <chr>  <dbl> <chr>   <chr>   <chr>    <chr>    <chr>    <chr>    <chr>    <chr>  <chr> <chr> <chr>         <dbl>
##  1 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…          72
##  2 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…         115
##  3 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…          96
##  4 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…          84
##  5 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…         114
##  6 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…         129
##  7 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…          36
##  8 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…          87
##  9 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…          74
## 10 PD63118  1     2488095 GCCT… G     TNFR…      1 .       .       .        .        .        1-10-de… .        no-SNV ENSP… del   plate1…         102
## # ℹ 151,994 more rows
## # ℹ 6 more variables: ref_depth <dbl>, mut_depth <dbl>, mut_vaf <dbl>, on_1p <lgl>, celltype_VDJ_SHM_CSR <fct>, loh_1p <chr>
```

We filter to `mut_vaf` > 0.25 and `mut_depth` > 1.


``` r
genos_filter <-
  genos %>%
  filter(mut_vaf > 0.2, mut_depth > 1)

# return
genos_filter
```

```
## # A tibble: 129 × 25
##    donor_id chr       pos ref   alt   gene  strand ref_cod mut_cod ref3_cod mut3_cod aachange ntchange codonsub impact pid   type  cell_id total_depth
##    <chr>    <chr>   <dbl> <chr> <chr> <chr>  <dbl> <chr>   <chr>   <chr>    <chr>    <chr>    <chr>    <chr>    <chr>  <chr> <chr> <chr>         <dbl>
##  1 PD63118  1     2488105 T     A     TNFR…      1 T       A       ATG      AAG      M1K      T2A      ATG>AAG  Misse… ENSP… snv   plate1…         131
##  2 PD63118  1     2488106 G     T     TNFR…      1 G       T       TGG      TTG      M1I      G3T      ATG>ATT  Misse… ENSP… snv   plate3…         147
##  3 PD63118  1     2488138 G     A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP… snv   plate1…         157
##  4 PD63118  1     2488138 G     A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP… snv   plate3…         109
##  5 PD63118  1     2488138 G     A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP… snv   plate3…          30
##  6 PD63118  1     2488138 G     A     TNFR…      1 G       A       TGG      TAG      W12*     G35A     TGG>TAG  Nonse… ENSP… snv   plate3…         134
##  7 PD63118  1     2488139 G     A     TNFR…      1 G       A       GGA      GAA      W12*     G36A     TGG>TGA  Nonse… ENSP… snv   plate1…          87
##  8 PD63118  1     2488139 G     A     TNFR…      1 G       A       GGA      GAA      W12*     G36A     TGG>TGA  Nonse… ENSP… snv   plate1…          94
##  9 PD63118  1     2489164 G     A     TNFR…      1 G       A       AGG      AAG      .        .        .        Essen… ENSP… snv   plate1…          31
## 10 PD63118  1     2489171 T     G     TNFR…      1 T       G       GTA      GGA      Y26D     T76G     TAT>GAT  Misse… ENSP… snv   plate1…          40
## # ℹ 119 more rows
## # ℹ 6 more variables: ref_depth <dbl>, mut_depth <dbl>, mut_vaf <dbl>, on_1p <lgl>, celltype_VDJ_SHM_CSR <fct>, loh_1p <chr>
```


``` r
n_muts <- genos_filter %>% distinct(chr, pos, ref, alt) %>% nrow()
n_muts_geno <- genos_filter %>% filter(mut_vaf > 0) %>% distinct(chr, pos, ref, alt) %>% nrow()
```

DNA and targeted DNA genotyping counts were piled up per cell. 40 / 1,434 
mutations were recapitulated in this population. All 40 mutations were found 
exclusively in activated B cells [see figure X]. The top mutated genes were 
TNFRSF14 (n mutations = 16, n cells mutated = 12; 6 missense, 5 nonsense, 3 
essential splice, 2 synonymous) and IGLL5 (n mutations = 16, n cells mutated = 
9; 12 missense, 3 synonymous, 1 essential splice). 

We plot the impact of the mutations genotyped.


``` r
genos_filter %>%
  distinct(chr, pos, ref, alt, impact, gene) %>%
  add_count(gene) %>%
  ggplot(aes(x = reorder(gene, -n), fill = impact)) +
  geom_bar() +
  geom_text(aes(label = n, y = n), vjust = -0.5) +
  scale_fill_manual(values = impact_colours) +
  guides(x = guide_axis(angle = -90)) +
  theme_classic() +
  labs(x = "", y = "n mutations")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_impact-1.png)<!-- -->

We generate a heatmap of all mutations.


``` r
plot_vaf_heatmap(genos, p_title = "all mutations", show_rownames = FALSE,
                 annotations = c("celltype_VDJ_SHM_CSR", "loh_1p"))
```


``` r
plot_vaf_heatmap(genos_filter, p_title = "filtered mutations",
                 show_rownames = FALSE,
                 annotations = c("celltype_VDJ_SHM_CSR", "loh_1p"))
```


``` r
p <-
  genos_filter %>%
  group_by(chr, pos, ref, alt) %>%
  filter(n() > 1) %>%
  plot_vaf_heatmap(p_title = "filtered shared mutations",
                   annotations = c("celltype_VDJ_SHM_CSR", "loh_1p"))
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_heatmap_filtered_shared-1.png)<!-- -->


``` r
p <-
  genos_filter %>%
  filter(impact != "Synonymous") %>%
  plot_vaf_heatmap("filtered non-synonymous mutations",
                   annotations = c("celltype_VDJ_SHM_CSR", "loh_1p"))
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_heatmap_filtered_nonsynonymous-1.png)<!-- -->


``` r
p <-
  genos_filter %>%
  filter(gene == "TNFRSF14") %>%
  plot_vaf_heatmap("non-synonymous TNFRSF14 mutations",
                   annotations = c("celltype_VDJ_SHM_CSR", "loh_1p"))
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_heatmap_filtered_TNFRSF14-1.png)<!-- -->


``` r
p <-
  genos_filter %>%
  filter(gene == "TNFRSF14", impact != "Synonymous") %>%
  plot_vaf_heatmap("TNFRSF14 mutations",
                   annotations = c("celltype_VDJ_SHM_CSR", "loh_1p"))
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_heatmap_filtered_TNFRSF14_nonsynonymous-1.png)<!-- -->

We plot TNFRSF14 mutant calls in the 1p LOH and non-1p LOH cells.


``` r
p <-
  genos_filter %>%
  filter(gene == "TNFRSF14") %>%
  ggplot(aes(x = loh_1p, y = mut_vaf)) +
  geom_boxplot(outlier.colour = NA) +
  geom_jitter(height = 0, width = 0.1) +
  theme_classic()
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

We plot the VAF distribution of the mutations.


``` r
genos %>%
  left_join(ss) %>%
  filter(mut_vaf > 0) %>%
  ggplot(aes(x = mut_vaf, fill = on_1p & loh_1p == "1")) +
  geom_histogram() +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

``` r
genos_filter %>%
  left_join(ss) %>%
  ggplot(aes(x = mut_vaf, fill = on_1p & loh_1p == "1")) +
  geom_histogram() +
  lims(x = c(0, NA)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette = "Set1")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/unnamed-chunk-5-2.png)<!-- -->


``` r
# % of B cells with a nonsynonymous mutation
n_mat_b_cells <-
  man_insp %>%
  filter(celltype_SHM == "mature B cell") %>%
  nrow()
p_dat <-
  genos_filter %>%
  filter(impact != "Synonymous") %>%
  group_by(gene) %>%
  summarise(n_cells = n_distinct(cell_id)) %>%
  mutate(prop_cells = n_cells / n_mat_b_cells)
p2 <-
  p_dat %>%
  ggplot(aes(x = reorder(gene, -prop_cells), y = prop_cells)) +
  geom_col(fill = "#d6630d") +
  theme_classic() +
  labs(x = "gene", y = "% of mature B cells with a driver mutation") +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = guide_axis(angle = -45))

# mutations per gene
p1 <-
  genos_filter %>%
  left_join(man_insp) %>%
  add_count(gene, name = "total") %>%
  count(gene, impact, total, celltype_SHM) %>%
  mutate(impact = ifelse(impact == "no-SNV", "Synonymous", impact) %>% factor(levels = rev(names(impact_colours)))) %>%
  left_join(p_dat %>% distinct(gene, prop_cells)) %>%
  ggplot(aes(x = reorder(gene, -prop_cells), y = n, fill = impact)) +
  geom_col() +
  scale_fill_manual(values = impact_colours) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(y = "n mutations")

# plot
p1 / p2
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_vaf_dist_hyb-1.png)<!-- -->


``` r
p_dat <-
  genos_filter %>%
  filter(impact != "Synonymous") %>%
  arrange(loh_1p) %>%
  mutate(cell_id = forcats::fct_inorder(cell_id)) %>%
  add_count(gene) %>%
  add_count(gene, cell_id, name = "n per cell") %>%
  mutate(gene = paste0(gene, " (", n, ")")) %>%
  arrange(n) %>%
  mutate(gene = forcats::fct_inorder(gene)) %>%
  distinct(cell_id, celltype_VDJ_SHM_CSR, gene, loh_1p, `n per cell`)

p1 <-
  p_dat %>%
  ggplot(aes(x = cell_id, y = 1, fill = loh_1p)) +
  geom_tile() +
  scale_fill_brewer(palette = "Set1") +
  theme_void()

p2 <-
  p_dat %>%
  ggplot(aes(x = cell_id, y = gene, fill = `n per cell`)) +
  geom_tile() +
  theme_classic() +
  guides(x = guide_axis(angle = -90)) +
  viridis::scale_fill_viridis() +
  coord_equal()

p1 / p2 + plot_layout(heights = c(1, 10), guides = "collect")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_drivers_per_cell-1.png)<!-- -->

# Somatic mutation calling

We load mutations generated by previous reports and we combine the calls.


``` r
# get nr and nv
bj_muts_38 <-
  c("dna", "dnahyb") %>%
  purrr::set_names() %>%
  purrr::map(function(seq_type) {
    c("mut_depth" = "NV", "total_depth" = "NR") %>%
      purrr::map(function(i) {
        paste0("out/BaseJumper/bj-somatic-variantcalling/", seq_type,
               "/PD63118_run/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA", 
               "/Sequoia_group_PD63118_bino-10_rhosnp0.4_rhoindel0.4_mincov10_maxcov500_both_",
               i, "_filtered_all.txt") %>%
          read.table() %>%
          tibble::rownames_to_column(var = "mut_id")
      }) %>%
      bind_rows(.id = "name") %>%
      pivot_longer(cols = -c("name", "mut_id"), names_to = "id") %>%
      pivot_wider() %>%
      separate_wider_delim(cols = "mut_id", delim = "_", cols_remove = FALSE,
                           names = c("chr", "pos", "ref", "mut")) %>%
      mutate(pos = as.numeric(pos)) %>%
      filter(mut_depth > 0)
  }) %>%
  bind_rows(.id = "seq_type") %>%
  left_join(ss %>% dplyr::select(id, cell_id)) %>%
  group_by(cell_id) %>%
  mutate(seq_type = paste(unique(seq_type), collapse = "+")) %>%
  # pileup across dna and dnahyb
  group_by(chr, pos, ref, mut, mut_id, cell_id, seq_type) %>%
  summarise(mut_depth = sum(mut_depth), total_depth = sum(total_depth)) %>%
  ungroup() %>%
  mutate(mut_vaf = mut_depth / total_depth)

# annotate bj_muts with dndscv
refcds <- "../reference/dndscv/RefCDS_human_GRCh38_GencodeV18_recommended.rda"
dndscv_in <-
  bj_muts_38 %>%
  transmute(sampleID = "", chr = gsub("chr", "", chr), pos, ref, mut) %>%
  distinct()
dndscv_out <-
  dndscv::dndscv(dndscv_in, max_muts_per_gene_per_sample = Inf,
                 max_coding_muts_per_sample = Inf, outp = 1, refdb = refcds)

# add mutation annotations
bj_muts_38 <-
  bj_muts_38 %>%
  left_join(
    dndscv_out$annotmuts %>%
      mutate(chr = paste0("chr", chr)) %>%
      dplyr::select(-sampleID)) %>%
  left_join(man_insp)

# lift over to grch37
bj_muts <-
  bj_muts_38 %>%
  {GRanges(seqnames = .$chr,
           ranges = IRanges(start = .$pos, end = .$pos),
           pos = .$pos)} %>%
  liftOver(import.chain("../reference/liftOver/hg38ToHg19.over.chain")) %>%
  unlist() %>%
  tibble::as_tibble() %>%
  distinct() %>%
  dplyr::transmute(chr = seqnames, pos_37 = start, pos) %>%
  inner_join(bj_muts_38) %>%
  mutate(pos = pos_37,
         mut_id = paste(chr, pos, ref, mut, sep = "_"))

# return
bj_muts
```

```
## # A tibble: 65,516 × 55
##    chr    pos_37     pos ref   mut   mut_id     cell_id seq_type mut_depth total_depth mut_vaf gene  strand ref_cod mut_cod ref3_cod mut3_cod aachange
##    <chr>   <int>   <int> <chr> <chr> <chr>      <chr>   <chr>        <int>       <int>   <dbl> <chr>  <dbl> <chr>   <chr>   <chr>    <chr>    <chr>   
##  1 chr1  1517661 1517661 G     A     chr1_1517… plate3… dna+dna…        19          33 0.576   <NA>      NA <NA>    <NA>    <NA>     <NA>     <NA>    
##  2 chr1  1761291 1761291 C     T     chr1_1761… plate3… dna+dna…         9          21 0.429   <NA>      NA <NA>    <NA>    <NA>     <NA>     <NA>    
##  3 chr1  1801468 1801468 A     AT    chr1_1801… plate3… dna+dna…         4           7 0.571   <NA>      NA <NA>    <NA>    <NA>     <NA>     <NA>    
##  4 chr1  2488094 2488094 A     C     chr1_2488… plate3… dna+dna…        53         114 0.465   <NA>      NA <NA>    <NA>    <NA>     <NA>     <NA>    
##  5 chr1  2488105 2488105 T     A     chr1_2488… plate1… dnahyb         151         151 1       TNFR…      1 T       A       ATG      AAG      M1K     
##  6 chr1  2488105 2488105 T     A     chr1_2488… plate1… dnahyb           1         100 0.01    TNFR…      1 T       A       ATG      AAG      M1K     
##  7 chr1  2488106 2488106 G     T     chr1_2488… plate3… dna+dna…        63         147 0.429   TNFR…      1 G       T       TGG      TTG      M1I     
##  8 chr1  2488138 2488138 G     A     chr1_2488… plate1… dnahyb           1         103 0.00971 TNFR…      1 G       A       TGG      TAG      W12*    
##  9 chr1  2488138 2488138 G     A     chr1_2488… plate1… dnahyb           2         191 0.0105  TNFR…      1 G       A       TGG      TAG      W12*    
## 10 chr1  2488138 2488138 G     A     chr1_2488… plate1… dnahyb          71         191 0.372   TNFR…      1 G       A       TGG      TAG      W12*    
## # ℹ 65,506 more rows
## # ℹ 37 more variables: ntchange <chr>, codonsub <chr>, impact <chr>, pid <chr>, n_cells <dbl>, plex <chr>, plate <dbl>, DNA_PrePCR_conc <dbl>,
## #   RNA_PrePCR_conc <dbl>, suspected_doublet <lgl>, doublet_rationale <chr>, celltype_bj_expression <chr>, celltype_SHM <chr>,
## #   celltype_VDJ_recomb <chr>, class_switch_recombination_CSR <chr>, chr_dropout <lgl>, loh_1p <chr>, TNFRSF14_mut <chr>, TNFRSF14_mut_VAF <dbl>,
## #   TNFRSF14_mut_in_NanoSeq_data <chr>, productive_heavy_chain <chr>, productive_heavy_chain_VDJ_consensus_sequences <chr>,
## #   productive_heavy_chain_CDR3_nt_IgBLAST <chr>, productive_heavy_chain_CDR3_nt_length <dbl>, productive_heavy_chain_CDR3_aa_IgBLAST <chr>,
## #   productive_heavy_chain_CDR3_aa_length <dbl>, unproductive_heavy_chain <chr>, unproductive_heavy_chain_VDJ_consensus_sequences <chr>, …
```

Plot the distribution of mutation types per cell.


``` r
# all muts
bj_muts %>%
  filter(!is.na(impact)) %>%
  mutate(group = "all") %>%
  bind_rows(bj_muts %>% filter(impact != "Non-coding") %>% mutate(group = "coding")) %>%
  add_count(group, cell_id) %>%
  ggplot(aes(x = reorder(cell_id, -n), fill = impact)) +
  geom_bar() +
  scale_fill_manual(values = impact_colours) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  guides(x = guide_axis(angle = -90)) +
  labs(x = "") +
  facet_grid(group ~ ., scales = "free_y")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250422_first_year_report_plots_files/figure-html/plot_mut_types_all-1.png)<!-- -->

Plot the VAF heatmap.


``` r
bj_muts %>%
  plot_vaf_heatmap(p_title = "all mutations", show_rownames = FALSE,
                   annotations = c("loh_1p"))
```

# VAF distribution

Plot the VAF distribution.


``` r
c(1, 5, 10, 20, 50) %>%
  purrr::walk(function(min_total_depth) {
    p_dat <- bj_muts %>% filter(total_depth >= min_total_depth)
    p <-
      p_dat %>%
      ggplot(aes(x = mut_vaf)) +
      geom_histogram() +
      theme_classic() +
      labs(title = paste("VAF distribution, min depth =", min_total_depth),
           subtitle = paste(nrow(p_dat), "mutations"), x = "alt VAF")
    print(p)
  })
```
