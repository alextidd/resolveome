---
title: "Shared mutations from the bj-somatic-variantcalling WGS output"
author: "Alexandra Tidd"
date: "21 March, 2025"
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
bj_dir <- "out/BaseJumper/bj-somatic-variantcalling/filter_cells/dnahyb/"
bj_run_dir <- file.path(bj_dir, "PD63118_250321_115728")
vcfs_subdir <- "/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SUBSET_VCF_VARIANTS/"
```

We load the samplesheet.


``` r
man_insp <-
  readr::read_tsv("data/manual_inspection/2024-12-20_PD63118_PTA_BAF_LoH_CellType_Mut_Summary.tsv") %>%
  dplyr::select(-id, -run_id)
metadata <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::left_join(man_insp) %>%
  dplyr::mutate(
    # define celltypes
    celltype = dplyr::case_when(
      celltype_VDJ_recomb == "alpha-beta T cell" ~ "T cell",
      celltype_VDJ_recomb == "B cell" & celltype_SHM == "mature B cell" ~
        "mature B cell",
      celltype_VDJ_recomb == "B cell" & celltype_SHM == "not mature B cell" ~
        "non-mature B cell",
      TRUE ~ "uncertain"),
    # define 1p_loh
    loh_1p = dplyr::case_when(`1p_loh` ~ "1p LOH", !`1p_loh` ~ "no 1p LOH",
                              TRUE ~ "unknown"))
ss <-
  readr::read_csv(file.path(bj_dir, "samplesheet.csv"))
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
  paste0(bj_run_dir, "/multiqc/multiqc_data/multiqc_PROJECT ID PLOT_TITLE_1.txt") %>%
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

# `BaseJumper` somatic variant calls

We load the variants and add cell-level metadata.


``` r
# get vafs, add cell-level metadata
vafs <- xfun::cache_rds({
  ss$biosampleName %>%
    purrr::set_names() %>%
    purrr::map(function(id) {
      vcf_file <-
        paste0(bj_run_dir, vcfs_subdir, id, "_somatic_filtered_variants.vcf.gz")
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
  labs(x = "")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_mut_types_all-1.png)<!-- -->


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

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_mut_types_coding-1.png)<!-- -->

# Mutation spectra

Get the mutational spectra in the cells.


``` r
# get trinucleotide context
vafs_trinuc <-
  xfun::cache_rds({
    vafs %>%
      get_mut_trinucs()
  }, file = "vafs_trinuc.rds", rerun = params$rerun)
```

## Plot trinucleotide spectrum (all cells)


``` r
plot_mut_trinucs(vafs_trinuc, "All cells")
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra-1.png)<!-- -->

## Plot trinucleotide spectrum per celltype


``` r
split(vafs_trinuc, vafs_trinuc$celltype) %>%
  purrr::walk(function(vafs_trinuc_i) {
    print(plot_mut_trinucs(vafs_trinuc_i, unique(vafs_trinuc_i$celltype)))
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_ct-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_ct-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_ct-3.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_ct-4.png)<!-- -->

## Plot trinucleotide spectrum per cell


``` r
split(vafs_trinuc, vafs_trinuc$id) %>%
  purrr::walk(function(vafs_trinuc_i) {
    print(plot_mut_trinucs(vafs_trinuc_i, unique(vafs_trinuc_i$id)))
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-3.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-4.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-5.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-6.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-7.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-8.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-9.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-10.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-11.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-12.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-13.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-14.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-15.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-16.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-17.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-18.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-19.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-20.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-21.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-22.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-23.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-24.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-25.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-26.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-27.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-28.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-29.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-30.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-31.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-32.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-33.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-34.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-35.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-36.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-37.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-38.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_spectra_per_cell-39.png)<!-- -->

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
           subtitle = paste(nrow(p_dat), "mutations"), x = "alt VAF")
    print(p)
  })
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_dist-1.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_dist-2.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_dist-3.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_dist-4.png)<!-- -->![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_dist-5.png)<!-- -->

# VAF heatmap

## Plot all mutations (all cells)


``` r
p <- plot_vaf_heatmap(vafs, "all muts")
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_heatmap-1.png)<!-- -->

## Plot all shared mutations (all cells)


``` r
p_dat_shared <-
  vafs %>%
  dplyr::mutate(n_cells = dplyr::n_distinct(id)) %>%
  dplyr::add_count(mut_id, name = "n_cells_w_mut") %>%
  dplyr::mutate(prop_cells_w_mut = n_cells_w_mut / n_cells) %>%
  # shared mutations
  dplyr::filter(n_cells_w_mut > 1)
p <- plot_vaf_heatmap(p_dat_shared, "all shared muts", show_rownames = TRUE)
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_heatmap_shared-1.png)<!-- -->

## Plot all shared mutations (% cells w mut < 0.3)

In order to exclude some of the germline noise, we look just at mutations 
present in <30% of cells.


``` r
p_dat_shared_filtered <-
  p_dat_shared %>%
  dplyr::filter(prop_cells_w_mut < 0.3)
p <- plot_vaf_heatmap(p_dat_shared_filtered,
                      "shared muts (% cells w mut < 0.3)", TRUE)
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_heatmap_shared_filtered-1.png)<!-- -->

## Plot coding shared mutations (all cells)


``` r
# look at coding muts
p_dat_coding <-
  p_dat_shared %>%
  dplyr::filter(impact != "Non-coding")
p <- plot_vaf_heatmap(p_dat_coding, "shared coding muts (all cells)", TRUE)
print(p)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_vaf_heatmap_shared_coding-1.png)<!-- -->

# Sequoia

We load the tree.


``` r
tree <-
  file.path(bj_run_dir, "/SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA/") %>%
  list.files(pattern = "_both_tree_with_branch_length.tree", full.names = TRUE) %>%
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

# colour celltypes
tip_annotations <-
  metadata %>%
  dplyr::filter(id %in% tree$tip.label) %>%
  dplyr::mutate(celltype = ifelse(is.na(celltype), "other", celltype)) %>% 
  {split(.$id, .$celltype)}
tip_annotations <- tip_annotations[tree$tip.label]

# plot tree with colored tip labels
plot(tree, tip.color = tip_col, cex = 0.6)
ape::axisPhylo(side = 1, backward = FALSE)
```

![](/lustre/scratch125/casm/team268im/at31/resolveome/reports/20250321_shared_muts_from_bj-somatic-variantcalling_dnahyb_filter_cells_files/figure-html/plot_tree-1.png)<!-- -->

# Driver mutations


``` r
vafs %>%
  dplyr::filter(gene %in% driver_genes) %>%
  dplyr::select(mut_id, gene, id, aachange, ntchange, codonsub, impact,
                ref_depth, alt_depth, alt_vaf) %>%
  knitr::kable()
```



|mut_id             |gene     |id                             |aachange |ntchange |codonsub |impact           | ref_depth| alt_depth|   alt_vaf|
|:------------------|:--------|:------------------------------|:--------|:--------|:--------|:----------------|---------:|---------:|---------:|
|chr1-2557736-T-C   |TNFRSF14 |plate3_wellA2_dnahyb_run49900  |L27P     |T80C     |CTC>CCC  |Missense         |         0|        44| 1.0000000|
|chr2-96144842-G-A  |DUSP2    |plate3_wellB2_dnahyb_run49900  |C143C    |C429T    |TGC>TGT  |Synonymous       |        26|        42| 0.6176471|
|chr1-2557792-G-T   |TNFRSF14 |plate3_wellB3_dnahyb_run49900  |E46*     |G136T    |GAG>TAG  |Nonsense         |         0|        12| 1.0000000|
|chr1-2556699-G-A   |TNFRSF14 |plate3_wellB4_dnahyb_run49900  |W12*     |G35A     |TGG>TAG  |Nonsense         |        67|        58| 0.4640000|
|chr1-2561748-C-G   |TNFRSF14 |plate3_wellB4_dnahyb_run49900  |S209R    |C627G    |AGC>AGG  |Missense         |        24|        36| 0.6000000|
|chr1-2557732-T-G   |TNFRSF14 |plate3_wellD4_dnahyb_run49900  |Y26D     |T76G     |TAT>GAT  |Missense         |         0|        20| 1.0000000|
|chr1-2557734-T-C   |TNFRSF14 |plate3_wellD4_dnahyb_run49900  |Y26Y     |T78C     |TAT>TAC  |Synonymous       |         0|        22| 1.0000000|
|chr1-2558453-C-T   |TNFRSF14 |plate3_wellD7_dnahyb_run49900  |Q97*     |C289T    |CAA>TAA  |Nonsense         |         0|        62| 1.0000000|
|chr1-2556699-G-A   |TNFRSF14 |plate3_wellD9_dnahyb_run49900  |W12*     |G35A     |TGG>TAG  |Nonsense         |         3|        25| 0.8928571|
|chr3-16377835-G-A  |RFTN1    |plate3_wellE10_dnahyb_run49900 |P237S    |C709T    |CCC>TCC  |Missense         |       266|       194| 0.4217391|
|chr6-31581809-C-G  |LTB      |plate3_wellE10_dnahyb_run49900 |.        |.        |.        |Essential_Splice |       153|        77| 0.3347826|
|chr1-2557732-T-G   |TNFRSF14 |plate3_wellE11_dnahyb_run49900 |Y26D     |T76G     |TAT>GAT  |Missense         |         0|        27| 1.0000000|
|chr1-2557734-T-C   |TNFRSF14 |plate3_wellE11_dnahyb_run49900 |Y26Y     |T78C     |TAT>TAC  |Synonymous       |         0|        27| 1.0000000|
|chr2-96144852-G-A  |DUSP2    |plate3_wellE11_dnahyb_run49900 |P140L    |C419T    |CCC>CTC  |Missense         |        68|        45| 0.3982301|
|chr9-5457168-T-G   |CD274    |plate3_wellE11_dnahyb_run49900 |L48V     |T142G    |TTA>GTA  |Missense         |       253|       332| 0.5675214|
|chr1-2556667-G-T   |TNFRSF14 |plate3_wellF11_dnahyb_run49900 |M1I      |G3T      |ATG>ATT  |Missense         |        81|        63| 0.4375000|
|chr1-2556699-G-A   |TNFRSF14 |plate3_wellF11_dnahyb_run49900 |W12*     |G35A     |TGG>TAG  |Nonsense         |        54|        82| 0.6029412|
|chr1-2559822-G-A   |TNFRSF14 |plate3_wellF6_dnahyb_run49900  |.        |.        |.        |Essential_Splice |         0|        34| 1.0000000|
|chr1-2558341-A-G   |TNFRSF14 |plate3_wellF7_dnahyb_run49900  |.        |.        |.        |Essential_Splice |         0|        64| 1.0000000|
|chr1-2557825-T-C   |TNFRSF14 |plate3_wellH5_dnahyb_run49900  |C57R     |T169C    |TGC>CGC  |Missense         |         0|        31| 1.0000000|
|chr3-183492153-T-C |KLHL6    |plate3_wellH5_dnahyb_run49900  |E547G    |A1640G   |GAG>GGG  |Missense         |        58|        81| 0.5827338|
|chr1-2561672-G-A   |TNFRSF14 |plate3_wellH6_dnahyb_run49900  |.        |.        |.        |Essential_Splice |         0|        69| 1.0000000|

# Beta binomial


``` r
# get nv and nr
mats <-
  list("NR", "NV") %>%
  purrr::set_names() %>%
  purrr::map(function(mat_type) {
    file.path(bj_run_dir, "SOMATIC_VARIANT_WORKFLOW_Heuristic_Filter_SEQUOIA") %>%
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
