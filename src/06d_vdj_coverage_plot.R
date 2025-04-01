#!/usr/bin/env Rscript
# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q yesterday -M80000 -R 'span[hosts=1] select[mem>80000] rusage[mem=80000]' -J vdj_cov -o log/%J_vdj_cov.out -e log/%J_vdj_cov.err 'module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=$HOME/R-tmp-4.4; Rscript src/06d_vdj_coverage_plot.R'

# libraries
library(magrittr)
library(ggplot2)
library(patchwork)
library(scales)
library(dplyr)
library(stringr)

# dirs
out_dir <- "out/vdj_coverage/filter_cells/"

# load samplesheet
ss <- readr::read_csv(file.path(out_dir, "samplesheet.csv"))

# load regions
ig_tcr_regions <-
  readr::read_tsv("out/vdj_coverage/regions/ig_tcr_regions.bed",
                  col_names = c("chr", "start", "end", "type")) %>%
  mutate(chr = as.character(chr), region = paste0(chr, "_", type))
ig_tcr_regions_pos <-
  ig_tcr_regions %>%
  group_by(chr, start, end, region, type) %>%
  reframe(pos = start:end)
ig_tcr_pos <-
  readr::read_tsv("out/vdj_coverage/regions/ig_tcr_genes.tsv") %>%
  mutate(chr = as.character(chr)) %>%
  group_by(chr, gene, start, end) %>%
  reframe(pos = start:end) %>%
  left_join(ig_tcr_regions_pos %>% select(chr, region, pos))
ig_tcr <-
  ig_tcr_pos %>%
  distinct(chr, gene, start, end, region) %>%
  mutate(
    segment = case_when(
      # Heavy chain (IGH)
      str_detect(gene, "^IGHV") ~ "V",
      str_detect(gene, "^IGHD") ~ "D",
      str_detect(gene, "^IGHJ") ~ "J",
      str_detect(gene, "^IGH[ACGEMD]") ~ "C",

      # Kappa light chain (IGK)
      str_detect(gene, "^IGKV") ~ "V",
      str_detect(gene, "^IGKJ") ~ "J",
      str_detect(gene, "^IGKC") ~ "C",

      # Lambda light chain (IGL)
      str_detect(gene, "^IGLV") ~ "V",
      str_detect(gene, "^IGLJ") ~ "J",
      str_detect(gene, "^IGLC") ~ "C",

      # TCR alpha (TRA)
      str_detect(gene, "^TRAV") ~ "V",
      str_detect(gene, "^TRAJ") ~ "J",
      str_detect(gene, "^TRAC") ~ "C",

      # TCR beta (TRB)
      str_detect(gene, "^TRBV") ~ "V",
      str_detect(gene, "^TRBJ") ~ "J",
      str_detect(gene, "^TRBC") ~ "C",

      # TCR delta (TRD)
      str_detect(gene, "^TRDV") ~ "V",
      str_detect(gene, "^TRDD") ~ "D",
      str_detect(gene, "^TRDJ") ~ "J",
      str_detect(gene, "^TRDC") ~ "C",

      # TCR gamma (TRG)
      str_detect(gene, "^TRGV") ~ "V",
      str_detect(gene, "^TRGJ") ~ "J",
      str_detect(gene, "^TRGC") ~ "C",

      TRUE ~ NA_character_  # Catch unexpected cases
    ) %>% factor(levels = c("V", "D", "J", "C")))

# make plots
purrr::walk2(ss$id, ss$donor_id, function(id_i, donor_id_i) {
  base_dir <- file.path(out_dir, donor_id_i, id_i)
  dir.create(file.path(base_dir, "plots"), showWarnings = FALSE)

  print(paste(donor_id_i, id_i))

  # check if the plots have already been create
  if(file.exists(paste0(base_dir, "/plots/", id_i, "_chr7_TCR_binned_cov.png"))) {
    print("plots already created")
  } else {
    print("load mean cov per gene")
    regions_file <- paste0(base_dir, "/mosdepth/", id_i, ".regions.bed.gz")
    dat <-
      readr::read_tsv(gzfile(regions_file),
                      col_names = c("chr", "start", "end", "gene", "mean_cov"),
                      show_col_types = FALSE) %>%
      mutate(id = id_i, chr = as.character(chr)) %>%
      left_join(ig_tcr) %>%
      mutate(gene = gene %>% forcats::fct_reorder(-start)) %>%
      {split(., .$region)}

    print("plot mean cov per gene")
    purrr::walk2(names(dat), dat, function(chr_i, chr_dat) {
      p <-
        chr_dat %>%
        ggplot(aes(x = gene, y = mean_cov, fill = segment)) +
        geom_col() +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        guides(x = guide_axis(angle = -90)) +
        theme_classic() +
        ggtitle(paste0(id_i, " - chr", chr_i, " genes - mean coverage")) +
        scale_fill_brewer(palette = "Dark2")
      ggsave(paste0(base_dir, "/plots/", id_i, "_chr", chr_i, "_mean_cov.png"),
             p, height = 5, width = 20)
    })

    # print("bin cov per 1kb")
    # per_base_file <- paste0(base_dir, "/mosdepth/", id_i,
    #                         ".per-base.ig_tcr_regions.bed.gz")
    # cov <-
    #   readr::read_tsv(gzfile(per_base_file),
    #                   col_names = c("chr", "start", "end", "cov", "chr2",
    #                                 "start_region", "end_region", "type"),
    #                   show_col_types = FALSE) %>%
    #   transmute(chr, start, end, cov, region = paste0(chr, "_", type),
    #             start_region, end_region) %>%
    #   # expand to all positions
    #   group_by(chr, start, end, region, start_region, end_region, cov) %>%
    #   reframe(pos = start:(end - 1)) %>%
    #   # just get positions within the region (remove overhangs from intersect)
    #   filter(pos >= start_region, pos < end_region) %>%
    #   # create 1kb bins
    #   mutate(bin = (pos - start_region) %/% 1000) %>%
    #   # bin per 1kb
    #   group_by(region, chr, bin, start_region, end_region) %>%
    #   summarise(mean_cov = mean(cov),
    #             start_bin = min(pos), end_bin = max(pos))
    # binned_file <- paste0(base_dir, "/mosdepth/", id_i,
    #                       ".1kb_binned_cov.ig_tcr_regions.tsv")
    # cov %>% readr::write_tsv(binned_file)

    # unique(cov$region) %>%
    #   purrr::set_names() %>%
    #   purrr::walk(function(region_i) {
    #     p1 <-
    #       cov %>%
    #       filter(region == region_i) %>%
    #       ggplot(aes(x = start_bin, y = mean_cov)) +
    #       geom_col(fill = "black", colour = "black") +
    #       theme_classic() +
    #       ggtitle(paste0(id_i, " - ", region_i, " - mean coverage per 1kb")) +
    #       labs(x = "") +
    #       scale_x_discrete(breaks = pretty_breaks(n = 20),
    #                       labels = scales::comma) +
    #       guides(x = guide_axis(angle = -90))
    #     p2 <-
    #       ig_tcr %>%
    #       filter(region == region_i) %>%
    #       mutate(mid = (start + end) / 2) %>%
    #       ggplot(aes(x = start, xend = end, y = 0)) +
    #       geom_segment(size = 3) +
    #       ggrepel::geom_text_repel(aes(x = mid, label = gene, y = 0), size = 2,
    #                               nudge_x = 0, max.overlaps = Inf, angle = -90) +
    #       theme_void()
    #     p <- p1 / p2 + plot_layout(heights = c(2, 1))

    #     # save
    #     ggsave(paste0(base_dir, "/plots/", id_i, "_chr", region_i,
    #                   "_binned_cov.png"),
    #            p, height = 5, width = 20)
    #  })
  }
})
