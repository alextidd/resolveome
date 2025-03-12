---
title: "Shared mutations from the bj-somatic-variantcalling WGS output"
author: "Alexandra Tidd"
date: "07 March, 2025"
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
  rerun: false
---



We pick out a BaseJumper run.


``` r
# get the latest run
bj_dir <- "out/BaseJumper/bj-somatic-variantcalling/wgs/PD63118_250306_232001/"
vcfs_subdir <- "/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SUBSET_VCF_VARIANTS/"
```

We load the samplesheet.


``` r
man_insp <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::mutate(run_id = as.character(run_id))
metadata <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::left_join(man_insp) %>%
  # define celltypes
  dplyr::mutate(
    celltype = dplyr::case_when(
      celltype_VDJ_recomb == "alpha-beta T cell" ~ "T cell",
      celltype_VDJ_recomb == "B cell" & celltype_SHM == "mature B cell" ~
        "mature B cell",
      celltype_VDJ_recomb == "B cell" & celltype_SHM == "not mature B cell" ~
        "non-mature B cell",
      TRUE ~ NA))
ss <-
  readr::read_csv("out/BaseJumper/bj-somatic-variantcalling/wgs/samplesheet.csv")
```

We load the driver genes.


``` r
driver_genes <-
  readLines("../trencadis-seq/data/thyroid/driver_genes/driver_genes.txt")
```

# `MultiQC`

We load the `multiqc` output.


``` r
multiqc <-
  paste0(bj_dir, "/multiqc/multiqc_data/multiqc_PROJECT ID PLOT_TITLE_1.txt") %>%
  readr::read_tsv() %>%
  dplyr::rename(id = Sample) %>%
  dplyr::left_join(metadata %>% dplyr::select(-bam)) 

multiqc %>%
  # clean
  dplyr::filter(!chr_dropout, !suspected_doublet) %>%
  ggplot(aes(x = reorder(id, -MEAN_COVERAGE), y = MEAN_COVERAGE)) +
  geom_col() +
  geom_errorbar(aes(ymin = MEAN_COVERAGE - SD_COVERAGE,
                    ymax = MEAN_COVERAGE + SD_COVERAGE),
                width = 0.2) +
  guides(x = guide_axis(angle = -90)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "mean coverage")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/multiqc-1.png)<!-- -->

``` r
multiqc %>%
  dplyr::left_join(metadata) %>%
  dplyr::mutate(chr_dropout = ifelse(chr_dropout, "chr dropout",
                                     "no chr dropout")) %>%
  ggplot(aes(x = reorder(id, -HET_SNP_SENSITIVITY), y = HET_SNP_SENSITIVITY)) +
  geom_col() +
  facet_grid(~ chr_dropout, scales = "free_x", space = "free_x") +
  guides(x = guide_axis(angle = -90)) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "heterozygous SNP sensitivity")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/multiqc-2.png)<!-- -->

# `BaseJumper` somatic variant calls

We load the variants and add cell-level metadata.


``` r
# get vafs, add cell-level metadata
vafs <- xfun::cache_rds({
  ss$biosampleName %>%
    purrr::set_names() %>%
    purrr::map(function(id) {
      vcf_file <-
        paste0(bj_dir, vcfs_subdir, id, "_somatic_filtered_variants.vcf.gz")
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
              tidyr::separate_longer_delim(cols = c("allele", "ad"),
                                           delim = ",") %>%
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
    dplyr::bind_rows(.id = "id") %>%
    # annotate celltypes
    dplyr::left_join(metadata)
  }, file = "vafs.rds", rerun = params$rerun)
```

# Variant annotation

We annotate the variants using `dndscv`.


``` r
refdb_file <-
  "../reference/dndscv/RefCDS_human_GRCh38_GencodeV18_recommended.rda"
vafs <-
  xfun::cache_rds({
    vafs %>%
      dplyr::transmute(sampleID = id, chr = gsub("chr", "", chr), pos, ref,
                       mut = alt) %>%
      dplyr::distinct() %>%
      dndscv(refdb = refdb_file, outp = 1,
             max_coding_muts_per_sample = Inf,
             max_muts_per_gene_per_sample = Inf) %>%
      {.$annotmuts} %>%
      dplyr::transmute(chr = paste0("chr", chr), pos, ref, alt = mut, gene, pid,
                       id = sampleID, aachange, ntchange, codonsub, impact) %>%
      dplyr::full_join(vafs) %>%
      dplyr::mutate(impact = ifelse(is.na(impact), "Non-coding", impact) %>%
                      factor(levels = rev(names(impact_colours))),
                    mut_id = paste(chr, pos, ref, alt, sep = "-"))
  }, file = "vafs_annotated.rds", rerun = params$rerun)
```

# Cell filtering

We mark cells with chromosomal dropout and suspected doublets.


``` r
vafs <-
  xfun::cache_rds({
    vafs %>%
      # mark cells with chr dropout and suspected doublets
      dplyr::mutate(qual = ifelse(!chr_dropout & !suspected_doublet,
                                  "good", "bad"))
  }, file = "vafs_w_filter.rds", rerun = params$rerun)
```

# Mutation types per cell

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
  labs(x = "") +
  facet_grid(~ qual, scales = "free_x", space = "free_x")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_mut_types_all-1.png)<!-- -->


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
  labs(x = "") +
  facet_grid(~ qual, scales = "free_x", space = "free_x")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_mut_types_coding-1.png)<!-- -->

# Mutation spectra

Get the mutational spectra in the good cells.


``` r
# get trinucleotide context
vafs_trinuc <-
  xfun::cache_rds({
    vafs %>%
      dplyr::filter(qual == "good") %>%
      get_mut_trinucs()
  }, file = "vafs_trinuc.rds", rerun = params$rerun)
```

## Plot trinucleotide spectrum (good cells)


``` r
plot_mut_trinucs(vafs_trinuc, "good cells")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_good-1.png)<!-- -->

## Plot trinucleotide spectrum per celltype (good cells)


``` r
split(vafs_trinuc, vafs_trinuc$celltype) %>%
  purrr::walk(function(vafs_trinuc_i) {
    print(plot_mut_trinucs(vafs_trinuc_i, unique(vafs_trinuc_i$celltype)))
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_good_per_ct-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_good_per_ct-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_good_per_ct-3.png)<!-- -->

## Plot trinucleotide spectrum per cell (good cells)


``` r
split(vafs_trinuc, vafs_trinuc$id) %>%
  purrr::walk(function(vafs_trinuc_i) {
    print(plot_mut_trinucs(vafs_trinuc_i, unique(vafs_trinuc_i$id)))
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-3.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-4.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-5.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-6.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-7.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-8.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-9.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-10.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-11.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-12.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-13.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-14.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-15.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-16.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-17.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-18.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-19.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-20.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-21.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-22.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-23.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-24.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-25.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-26.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-27.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-28.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-29.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-30.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-31.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-32.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-33.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-34.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-35.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-36.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-37.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-38.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-39.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-40.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-41.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-42.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_per_cell-43.png)<!-- -->

# VAF distribution

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
           subtitle = paste(nrow(p_dat), "mutations"), x = "alt VAF") +
      facet_grid(qual ~ .)
    print(p)
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_dist-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_dist-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_dist-3.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_dist-4.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_dist-5.png)<!-- -->

# VAF heatmap

## Plot all mutations (all cells)


``` r
p <- plot_vaf_heatmap(vafs, "all muts")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_heatmap-1.png)<!-- -->

## Plot all shared mutations (all cells)


``` r
p_dat_shared <-
  vafs %>%
  dplyr::mutate(n_cells = dplyr::n_distinct(id)) %>%
  dplyr::add_count(mut_id, name = "n_cells_w_mut") %>%
  dplyr::mutate(prop_cells_w_mut = n_cells_w_mut / n_cells) %>%
  # shared mutations
  dplyr::filter(n_cells_w_mut > 1)
p <- plot_vaf_heatmap(p_dat_shared, "all shared muts")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_heatmap_shared-1.png)<!-- -->

## Plot all shared mutations (good cells)


``` r
p_dat_shared_good <- p_dat_shared %>% dplyr::filter(qual == "good")
p <- plot_vaf_heatmap(p_dat_shared_good, "all shared muts (good cells)")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_heatmap_shared_good-1.png)<!-- -->

## Plot all shared mutations (good cells, % cells w mut < 0.3)

In order to exclude some of the germline noise, we look just at mutations 
present in <30% of cells.


``` r
p_dat_shared_filtered <-
  p_dat_shared_good %>%
  dplyr::filter(prop_cells_w_mut < 0.3)
p <- plot_vaf_heatmap(p_dat_shared_filtered,
                      "shared muts (good cells, % cells w mut < 0.3)")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_heatmap_shared_good_filtered-1.png)<!-- -->

## Plot coding shared mutations (all cells)


``` r
# look at coding muts
p_dat_coding <-
  p_dat_shared %>%
  dplyr::filter(impact != "Non-coding")
p <- plot_vaf_heatmap(p_dat_coding, "shared coding muts (all cells)", TRUE)
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_vaf_heatmap_shared_coding-1.png)<!-- -->

# Case study: `plate3_wellD4` & `plate3_wellE11`

We see that cells `plate3_wellD4` and `plate3_wellE11` share a block of 
mutations. We will look at the spectra of these mutations that are exclusively 
shared between these two cells.


``` r
excl_vafs_trinuc <-
  vafs %>%
  dplyr::group_by(chr, pos, ref, alt) %>%
  dplyr::filter(paste(sort(unique(id)), collapse = ",") == "plate3_wellD4_dna_run49882,plate3_wellE11_dna_run49882") %>%
  get_mut_trinucs()
plot_mut_trinucs(excl_vafs_trinuc,
                 "plate3_wellD4, plate3_wellE11 exclusive shared muts")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_spectra_shared_exclusive-1.png)<!-- -->

# Sequoia

We load the tree.


``` r
tree <-
  file.path(bj_dir, "/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA/") %>%
  list.files(pattern = "tree_with_branch_length.tree", full.names = TRUE) %>%
  ape::read.tree()
```

## Phylogeny (all cells)


``` r
# colour clonal tipes
clones <-
  list(
    clone1 = c("plate3_wellD4_dna_run49882", "plate3_wellE11_dna_run49882"),
    clone2 = c("plate3_wellD11_dna_run49882", "plate3_wellA10_dna_run49882"),
    clone3 = c("plate3_wellD3_dna_run49882", "plate3_wellB7_dna_run49882"),
    clone4 = c("plate3_wellE8_dna_run49882", "plate3_wellF10_dna_run49882"))
tip_col <- rep("black", length(tree$tip.label)) %>% setNames(tree$tip.label)
tip_pal <- RColorBrewer::brewer.pal(length(clones), "Dark2")
purrr::walk2(clones, tip_pal, function(clone, col) {
  tip_col[clone] <<- col
})

# plot tree with colored tip labels
plot(tree, tip.color = tip_col, cex = 0.6)
ape::axisPhylo(side = 1, backward = FALSE)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_tree-1.png)<!-- -->

## Phylogeny (good cells)


``` r
good_ids <- dplyr::filter(vafs, qual == "good")$id %>% unique()
good_tree <- ape::drop.tip(tree, setdiff(tree$tip.label, good_ids))
tip_col <-
  rep("black", length(good_tree$tip.label)) %>%
  setNames(good_tree$tip.label)
purrr::walk2(clones, tip_pal, function(clone, col) {
  tip_col[clone] <<- col
})

plot(good_tree, tip.color = tip_col, cex = 0.6)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250304_plot_shared_muts_from_bj-somatic-variantcalling_wgs_files/figure-html/plot_tree_good-1.png)<!-- -->

# Driver mutations


``` r
vafs %>%
  dplyr::filter(gene %in% driver_genes) %>%
  dplyr::select(mut_id, gene, id, aachange, ntchange, codonsub, impact,
                ref_depth, alt_depth, alt_vaf) %>%
  knitr::kable()
```



|mut_id             |gene  |id                          |aachange |ntchange |codonsub |impact   | ref_depth| alt_depth|   alt_vaf|
|:------------------|:-----|:---------------------------|:--------|:--------|:--------|:--------|---------:|---------:|---------:|
|chr3-16377835-G-A  |RFTN1 |plate3_wellE10_dna_run49882 |P237S    |C709T    |CCC>TCC  |Missense |        18|        19| 0.5135135|
|chr9-5457255-C-T   |CD274 |plate3_wellE5_dna_run49882  |Q77*     |C229T    |CAG>TAG  |Nonsense |         0|        19| 1.0000000|
|chr6-167137290-G-T |CCR6  |plate3_wellE8_dna_run49882  |E354*    |G1060T   |GAA>TAA  |Nonsense |        12|        18| 0.6000000|
|chr6-167137290-G-T |CCR6  |plate3_wellF10_dna_run49882 |E354*    |G1060T   |GAA>TAA  |Nonsense |         0|        25| 1.0000000|
|chr6-31581846-G-A  |LTB   |plate3_wellH5_dna_run49882  |A59V     |C176T    |GCC>GTC  |Missense |         5|         6| 0.5454545|

# Beta binomial


``` r
# get nv and nr
mats <-
  list("NR", "NV") %>%
  purrr::set_names() %>%
  purrr::map(function(mat_type) {
    file.path(bj_dir, "SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA") %>%
      list.files(pattern = paste0("Mat_", mat_type), full.names = TRUE) %>%
      readr::read_tsv(col_names = FALSE) %>%
      tibble::column_to_rownames("X1") %>%
      as.matrix()
  })

# estimate overdispersion
estimateRho_gridml <- function(NV_vec, NR_vec) {
  # rho will be bounded within 1e-6 and 0.89
 rhovec <- 10^seq(-6, -0.05, by = 0.05)
 mu <- sum(NV_vec) / sum(NR_vec)
 ll <- sapply(rhovec, function(rhoj) sum(dbetabinom(x = NV_vec, size = NR_vec, rho = rhoj, prob = mu, log = TRUE)))
 return(rhovec[ll == max(ll)][1])
}

rho_est <- pvalue <- rep(NA, nrow(mats$NR))
require(VGAM)
for (k in 1:nrow(mats$NR)) {
	rho_est[k] <- estimateRho_gridml(NV_vec = as.numeric(mats$NV[k,]), NR_vec = as.numeric(mats$NR[k,]))
	if (k%%1000 == 0) {
		print(k)
	}
}

# num_samples: those with at least 2 supporting reads
df <-
  data.frame(pos = 1:length(rho_est),
             rho_est = log(rho_est),
             num_samples = rowSums(NV[rownames(NV), ] > 1))

ggplot(df, aes(x = pos, y = rho_est)) +
  geom_point(aes(colour = factor(num_samples), alpha = 0.5)) +
  theme(legend.position = "none")

hist(df$num_samples, breaks = 20, xlab = "Number of samples with mutation",
     ylab = "Number of mutations", main = "pmin(20)")

df$global_vaf <-
  apply(alts[rownames(vafs),], 1, sum) / apply(covs[rownames(vafs), ], 1, sum)
ggplot(df, aes(x = pos, y = rho_est)) +
  geom_point(aes(colour = global_vaf))
```
