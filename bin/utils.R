# function: plot VAF heatmap
plot_vaf_heatmap <- function(p_dat, p_title = "", annotations = c(), show_rownames = TRUE) {
  library(magrittr)

  # prepare data
  p_dat2 <-
    p_dat %>%
    dplyr::ungroup() %>%
    # only keep mutations with VAF > 0
    dplyr::filter(mut_vaf > 0) %>%
    dplyr::mutate(mut_id = paste(gene, chr, pos, ref, alt, sep = "_")) %>%
    # count number of cells with each mutation
    dplyr::add_count(mut_id, name = "n_cells_w_mut")

  # reshape data for heatmap
  heatmap_data <-
    p_dat2 %>%
    reshape2::dcast(mut_id + n_cells_w_mut ~ cell_id, value.var = "mut_vaf")
  rownames(heatmap_data) <- heatmap_data$mut_id
  heatmap_matrix <-
    heatmap_data %>%
    dplyr::select(-mut_id, -n_cells_w_mut) %>%
    as.matrix()
  n_total_muts <- dplyr::n_distinct(p_dat2$mut_id)
  n_total_cells <- dplyr::n_distinct(p_dat2$cell_id)

  # calculate the number of mutations per column (id) + other annots
  annotations_col <-
    p_dat2 %>%
    dplyr::count(cell_id, !!!syms(annotations)) %>%
    tibble::column_to_rownames("cell_id")

  # replace NAs with 0
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # prepare annotations for columns
  annotation_row <- data.frame(n_cells_w_mut = heatmap_data$n_cells_w_mut)
  rownames(annotation_row) <- heatmap_data$mut_id

  # plot
  pheatmap::pheatmap(
    mat = heatmap_matrix,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = show_rownames,
    annotation_row = annotation_row,
    annotation_col = annotations_col,
    color = rev(viridis::magma(100)),
    main = paste(p_title, "\nVAF heatmap\n",
                 n_total_cells, "cells,", n_total_muts, "mutations")
  )
}

# function: plot vaf distribution with depth information
plot_vaf_dist <- function(p_dat, p_title = "VAF distribution") {
  n_muts <-
    p_dat %>%
    dplyr::distinct(chr, pos, ref, mut) %>%
    nrow()
  p_dat_binned <-
    p_dat %>%
    dplyr::mutate(
      mut_vaf = ifelse(total_depth == 0, 0, mut_vaf),
      total_depth_bin = cut(
        total_depth, breaks = c(-Inf, 0, 10, 20, 50, 100, Inf),
        labels = c("0", "<10", "10–20", "20–50", "50–100", ">100")),
      mut_vaf_bin = cut(mut_vaf, breaks = seq(0, 1, length.out = 30 + 1),
                        include.lowest = TRUE),
      mut_vaf_min = stringr::str_extract(mut_vaf_bin, "[\\d\\.]+") %>% as.numeric(),
      mut_vaf_high = as.numeric(stringr::str_extract(mut_vaf_bin, "(?<=,)\\s*[\\d\\.]+")),
      mut_vaf_mid = (mut_vaf_min + mut_vaf_high) / 2) %>%
    dplyr::count(mut_vaf_bin, mut_vaf_mid, total_depth_bin)
  p <-
    p_dat_binned %>%
    ggplot2::ggplot(ggplot2::aes(x = mut_vaf_mid, y = n, fill = total_depth_bin)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::theme_classic() +
    ggplot2::labs(title = p_title, subtitle = paste(n_muts, "mutations"),
         x = "VAF bin", y = "n mutations") +
    ggplot2::lims(x = c(0, 1)) +
    ggplot2::scale_y_continuous(expand = c(0, 0))
  print(p)
}
