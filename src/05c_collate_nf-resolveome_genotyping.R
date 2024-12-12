# libraries
library(magrittr)

# load ss
ss <- 
  system("ls out/nf-resolveome/*/samplesheet.tsv", intern = TRUE) %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows() %>%
  dplyr::filter(run %in% c("49882", "49900"))

# load mutations
muts <-
  readr::read_tsv("out/nf-resolveome/annotated_mutations.tsv") %>%
  dplyr::filter(source == "nanoseq_mutations")

# remove dodgy sites identified across all Exome NanoSeq data via Andrew's Kolmogorov-Smirnoff test
bad_sites <- 
  "data/nanoseq/2024-11-13_ExomeQposRecurrentMutsKolmogorovSmirnov_Indels.tsv" %>%
  data.table::fread(data.table = FALSE, header = TRUE) %>%
  dplyr::filter(ks_test_q < 0.05)

# load genotyping
geno <-
  ss %>%
  purrr::pmap(function(run, id, donor_id, ...) {
    paste0("out/nf-resolveome/", run, "/", donor_id, "/", id, "/genotyping/",
           id, "_genotyped_mutations.tsv") %>%
      readr::read_tsv() %>%
      dplyr::inner_join(muts) %>%
      dplyr::mutate(run = run, id = id, donor_id = donor_id)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(mut_vaf = mut_depth / total_depth)

# save 
pos_geno <-
  geno %>%
  dplyr::filter(mut_depth > 5, mut_vaf > 0.1) %>%
  dplyr::select(mut_vaf, everything()) %>%
  dplyr::arrange(-mut_vaf, -mut_depth) %>%
  dplyr::mutate(mut_id = paste(chr, pos, ref, mut, sep = "_"),
                mut_lab = paste(chr, pos, ref, mut, type, aachange, ifelse(is.na(gene), "", gene))) %>%
  # remove bad sites
  dplyr::filter(!mut_id %in% bad_sites$Var1)
pos_geno %>%
  readr::write_tsv("out/nf-resolveome/annotated_genotyped_mutations.tsv")

# Reshape the data to a matrix format with mut_lab as rows and id as columns
plot_vaf_heatmap <- function(p_dat, p_title = "Mutational VAF Heatmap") {
  heatmap_data <-
    p_dat %>%
    dplyr::distinct() %>%
    dplyr::group_by(mut_lab) %>%
    dplyr::mutate(n_cells_w_mut = dplyr::n_distinct(id)) %>%
    dplyr::ungroup() %>%
    reshape2::dcast(mut_lab + n_cells_w_mut ~ id, value.var = "mut_vaf")
  rownames(heatmap_data) <- heatmap_data$mut_lab
  heatmap_matrix <- as.matrix(heatmap_data[ , -c(1:3)]) # Drop annotations

  # Calculate the number of mutations per column (id)
  n_muts <-
    p_dat %>%
    dplyr::group_by(id) %>%
    dplyr::summarise(n_muts = dplyr::n_distinct(mut_lab))
  annotation_col <- data.frame(n_muts = n_muts$n_muts)
  rownames(annotation_col) <- n_muts$id

  # Replace NAs with 0 or another value depending on your context
  heatmap_matrix[is.na(heatmap_matrix)] <- 0

  # Prepare annotations for columns
  annotation_row <- data.frame(n_cells_w_mut = heatmap_data$n_cells_w_mut)
  rownames(annotation_row) <- heatmap_data$mut_lab

  pheatmap::pheatmap(
    mat = heatmap_matrix,
    cluster_rows = TRUE,        # Perform hierarchical clustering on rows
    cluster_cols = TRUE,        # Perform hierarchical clustering on columns
    annotation_row = annotation_row, # Add row annotations
    annotation_col = annotation_col, # Add column annotations
    color = colorRampPalette(c("white", "blue"))(50), # Customize the color gradient
    main = p_title
  )
}

# Create the heatmap with row annotations
pdf("test.pdf", width = 20, height = 20)
pos_geno %>% dplyr::filter(run == "49900") %>% plot_vaf_heatmap()
pos_geno %>% dplyr::filter(run == "49900", !is.na(gene)) %>% plot_vaf_heatmap("Mutational VAF Heatmap (dndscv-annotated genes only)")
dev.off()