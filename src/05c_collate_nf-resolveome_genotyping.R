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
  dplyr::arrange(-mut_vaf, -mut_depth)
pos_geno %>%
  readr::write_tsv("out/nf-resolveome/annotated_genotyped_mutations.tsv")

geno %>%
  dplyr::filter(mut_depth > 1) %>%
  dplyr::select(mut_vaf, everything()) %>%
  dplyr::arrange(-mut_vaf, -mut_depth) 
