# libraries
library(magrittr)
library(biomaRt)
library(pheatmap)
library(ggplot2)

# get featurecounts
fc <-
  "out/featureCounts/PD63118_49901.featureCounts.txt" %>%
  readr::read_tsv(comment = "#")

# get summary
fc_summ <-
  "out/featureCounts/PD63118_49901.featureCounts.txt.summary" %>%
  readr::read_tsv() %>%
  dplyr::rename_with(~ gsub("\\..*", "", basename(.x)),
                     dplyr::starts_with("data/")) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("plex"),
                      names_to = "id") %>%
  dplyr::mutate(Status = factor(Status, levels = rev(sort(unique(Status))))) %>%
  dplyr::group_by(id) %>%
  dplyr::mutate(total = sum(value))

# plot summary
p1 <-
  fc_summ %>%
  ggplot(aes(x = reorder(id, total), y = value, fill = Status)) +
  geom_col(position = "stack") +
  coord_flip() +
  labs(x = "", y = "n reads") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.vector(pals::alphabet())) +
  ggtitle("featureCounts summary, n reads")
p2 <-
  fc_summ %>%
  ggplot(aes(x = reorder(id, total), y = value, fill = Status)) +
  geom_col(position = "fill") +
  coord_flip() +
  labs(x = "", y = "prop reads") +
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = as.vector(pals::alphabet())) +
  ggtitle("featureCounts summary, prop reads")

# get celltype curated markers from celltypist
ct_markers <-
  readr::read_tsv("../reference/celltypist/Basic_celltype_information.tsv") %>%
  janitor::clean_names() %>%
  tidyr::separate_longer_delim("curated_markers", delim = ", ")

# convert ensembl gene ids to gene symbols
ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
                      version = 104)
results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = gsub("\\..*", "", fc$Geneid),
                 mart = ensembl)

# add gene symbols
fc_annot <-
  fc %>%
  dplyr::mutate(ensembl_gene_id = gsub("\\..*", "", Geneid)) %>%
  dplyr::left_join(results) %>%
  dplyr::mutate(hgnc_symbol = ifelse(hgnc_symbol == "", NA, hgnc_symbol)) %>%
  dplyr::rename_with(~ gsub("\\..*", "", basename(.x)),
                     dplyr::starts_with("data/")) %>%
  # collapse duplicate rows with the same counts
  dplyr::select(hgnc_symbol, dplyr::starts_with("plex")) %>%
  dplyr::distinct() %>%
  # sum duplicate gene rows
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::summarise(dplyr::across(dplyr::starts_with("plex"), sum))

# add celltypes
fc_annot_w_cts <-
  fc_annot %>%
  dplyr::inner_join(
    ct_markers %>%
      dplyr::distinct(ct = high_hierarchy_cell_types,
                      hgnc_symbol = curated_markers)) %>%
  dplyr::arrange(ct) %>%
  dplyr::mutate(gene = make.unique(hgnc_symbol))

# create matrix for heatmap
fc_mat <-
  fc_annot_w_cts %>%
  dplyr::select(gene, dplyr::starts_with("plex")) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# normalise matrix
fc_mat_norm <- log2(fc_mat + 1)

# row annotations
row_annotation <-
  fc_annot_w_cts %>%
  dplyr::select(gene, ct) %>%
  tibble::column_to_rownames("gene")

# heatmap normalised
p3 <-
  pheatmap(
  fc_mat_norm,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  annotation_row = row_annotation,
  labels_row = fc_annot_w_cts$hgnc_symbol,
  fontsize = 10, fontsize_row = 7, fontsize_col = 8,
  main = "log2(count + 1) heatmap"
)

# heatmap unnormalised
p4 <-
  pheatmap(
  fc_mat,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  annotation_row = row_annotation,
  labels_row = fc_annot_w_cts$hgnc_symbol,
  fontsize = 10, fontsize_row = 7, fontsize_col = 8,
  main = "count heatmap"
)

# save
pdf("out/featureCounts/PD63118_49901.featureCounts.summary.pdf",
    width = 12, height = 10)
print(p1)
print(p2)
dev.off()

pdf("out/featureCounts/PD63118_49901.featureCounts.celltype_heatmap.pdf",
    width = 15, height = 20)
print(p3)
dev.off()


pdf("out/featureCounts/PD63118_49901.featureCounts.celltype_heatmap_counts.pdf",
    width = 15, height = 20)
print(p4)
dev.off()