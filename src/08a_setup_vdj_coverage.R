# libraries
library(magrittr)
library(ggplot2)

# load data
df <-
  readr::read_tsv("data/vdj_reconstruction/ig_tcr_genes_pseudogenes.tsv") %>%
  janitor::clean_names() %>%
  dplyr::select(chr = chromosome_scaffold_name, start = gene_start_bp,
                end = gene_end_bp, gene = gene_name, everything()) %>%
  dplyr::mutate(class = substr(gene, 1, 4)) %>%
  dplyr::filter(chr %in% c("2", "7", "14", "22")) 

# write tsv
df %>%
  readr::write_tsv("out/vdj_reconstruction/ig_tcr_genes.tsv")

# write bed
df %>%
  dplyr::select(chr, start, end, gene) %>%
  readr::write_tsv("out/vdj_reconstruction/ig_tcr_genes.bed", col_names = FALSE)

# write samplesheet
readr::read_csv("data/resolveome/DNA/samplesheet_local.csv") %>%
  dplyr::filter(run == 49882) %>%
  readr::write_csv("out/vdj_reconstruction/samplesheet.csv")
