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
    segment = substr(gene, 4, 4),
    segment = ifelse(segment %in% c("A", "E", "G", "M"), "C", segment) %>%
      factor(levels = c("V", "D", "J", "C")))

# make plots
purrr::walk2(ss$id, ss$donor_id, function(id_i, donor_id_i) {
  base_dir <- file.path(out_dir, donor_id_i, id_i)
  dir.create(file.path(base_dir, "plots"), showWarnings = FALSE)

  print(paste(donor_id_i, id_i))

  print("load mean cov per gene")
  regions_file <- paste0(base_dir, "/mosdepth/", id_i, ".regions.bed.gz")
  dat <-
    readr::read_tsv(gzfile(regions_file),
                    col_names = c("chr", "start", "end", "gene", "mean_cov"),
                    show_col_types = FALSE) %>%
    mutate(id = id_i, chr = as.character(chr)) %>%
    left_join(ig_tcr) %>%
    mutate(gene = gene %>% forcats::fct_reorder(start)) %>%
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
})
