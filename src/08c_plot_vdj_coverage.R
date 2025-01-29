# libraries
library(magrittr)
library(ggplot2)

# load samplesheet
ss <-
  readr::read_csv("out/vdj_reconstruction/samplesheet.csv")

# load regions
ig_tcr <-
  readr::read_tsv("out/vdj_reconstruction/regions/ig_tcr_genes.tsv") %>%
  dplyr::mutate(chr = as.character(chr))
ig_tcr_pos <-
  ig_tcr %>%
  dplyr::group_by(chr, gene, start, end) %>%
  dplyr::reframe(pos = start:end)

# make plots
purrr::walk2(ss$id, ss$donor_id, function(id_i, donor_id) {
  base_dir <- paste0("out/vdj_reconstruction/", donor_id, "/", id_i)
  dir.create(paste0(base_dir, "/plots"), showWarnings = FALSE)

  print(paste(donor_id, id_i))

  print("load mean cov per gene")
  regions_file <- paste0(base_dir, "/mosdepth/", id_i, ".regions.bed.gz")
  dat <-
    readr::read_tsv(gzfile(regions_file),
                    col_names = c("chr", "start", "end", "gene", "mean_cov"),
                    show_col_types = FALSE) %>%
    dplyr::mutate(id = id_i, chr = as.character(chr),
                  class = dplyr::case_when(grepl("^IG", gene) ~ "BCR",
                                           grepl("^TR", gene) ~ "TCR",
                                           TRUE ~ NA),
                  chr_x_class = paste0(chr, "_", class)) %>%
    dplyr::left_join(ig_tcr) %>%
    dplyr::mutate(gene = gene %>% forcats::fct_reorder(start)) %>%
    {split(., .$chr_x_class)}

  print("plot mean cov per gene")
  purrr::walk2(names(dat), dat, function(chr_i, chr_dat) {
    p <-
      chr_dat %>%
      ggplot(aes(x = gene, y = mean_cov)) +
      geom_col(fill = "black") +
      scale_x_discrete(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      guides(x = guide_axis(angle = -90)) +
      theme_classic() +
      ggtitle(paste0(id_i, " - chr", chr_i, " genes - mean coverage"))
    ggsave(paste0(base_dir, "/plots/", id_i, "_chr", chr_i, "_mean_cov.png"),
           p, height = 5, width = 20)
  })

  print("load cov per base")
  per_base_file <- paste0(base_dir, "/mosdepth/", id_i,
                          ".per-base.regions.bed.gz")
  cov <-
    readr::read_tsv(gzfile(per_base_file),
                    col_names = c("chr", "start", "end", "cov", "chr2",
                                  "start_gene", "end_gene", "gene"),
                    show_col_types = FALSE) %>%
    dplyr::mutate(id = id_i, chr = as.character(chr),
                  gene = gene %>% forcats::fct_reorder(start),
                  class = dplyr::case_when(grepl("^IG", gene) ~ "BCR",
                                           grepl("^TR", gene) ~ "TCR",
                                           TRUE ~ NA),
                  chr_x_class = paste0(chr, "_", class)) %>%
    {split(., .$chr_x_class)}

  print("plot cov per base")
  purrr::walk2(names(cov), cov, function(chr_i, chr_dat) {
    overall_min <- min(chr_dat$start)
    overall_max <- max(chr_dat$end)
    p <-
      chr_dat %>%
      ggplot() +
      geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = cov),
                fill = "black") +
      scale_x_continuous(breaks = c(overall_min, overall_max),
                         expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      guides(x = guide_axis(angle = 90)) +
      facet_grid(~ gene, scale = "free_x", space = "free_x") +
      theme_classic() +
      theme(strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
            panel.spacing = unit(0, "lines"),
            panel.border = element_rect(color = "grey", fill = NA,
                                        linewidth = 0.1),
            strip.background = element_rect(color = "grey", fill = NA,
                                            linewidth = 0.1,
                                            linetype = "solid")) +
      ggtitle(paste0(id_i, " - chr", chr_i, " genes - coverage per base"))
    ggsave(paste0(base_dir, "/plots/", id_i, "_chr", chr_i, "_cov.png"),
           p, height = 5, width = 20)
  })
})
