# libraries
library(magrittr)

# dirs
donor_id <- "PD63118"
out_dir <- paste0("out/nf-resolveome/", donor_id, "/")
dir.create(out_dir)
wd <- getwd()

# genotype muts and snps
muts_and_snps <-
  c("caveman_snps", "nanoseq_mutations") %>%
  purrr::set_names() %>%
  purrr::map(function(source) {
    paste0("out/genotyping/", source, ".tsv") %>%
      readr::read_tsv()
  }) %>%
  dplyr::bind_rows(.id = "source")
muts_and_snps %>%
  readr::write_tsv(paste0(out_dir, "/mutations.tsv"))

# generate samplesheets
ss <- list()

# samplesheet 49686 (19 cells)
ss[["49686"]] <-
  tibble::tibble(
    run = "49686",
    lane_n = "4-5",
    plex_n = seq(1, 19)) %>%
  dplyr::transmute(
    run, plex_n, lane_n,
    id = paste0(run, "_plex", plex_n),
    bam = paste0("/seq/illumina/runs/", substr(run, 1, 2), "/", run, "/",
                 "lane", lane_n, "/plex", plex_n, "/", run, "_", lane_n, "#",
                 plex_n, ".cram"))

# samplesheet 49882 (80 cells)
ss[["49882"]] <-
  tibble::tibble(
    run = "49882",
    plex_n = seq(1, 80)) %>%
  dplyr::transmute(
    run, plex_n,
    id = paste0(run, "_plex", plex_n),
    bam = paste0("/seq/illumina/runs/", substr(run, 1, 2), "/", run, "/plex",
                 plex_n, "/", run, "#", plex_n, ".cram"))

# samplesheet 49900 (80 cells with bait capture)
ss[["49900"]] <-
  tibble::tibble(
    run = "49900", 
    plex_n = seq(1, 80),
    lane_n = "2") %>%
  dplyr::transmute(
    run, plex_n, lane_n,
    id = paste0(run, "_plex", plex_n),
    bam = paste0("/seq/illumina/runs/", substr(run, 1, 2), "/", run, "/lane",
                 lane_n, "/plex", plex_n, "/", run, "_", lane_n, "#", plex_n,
                 ".cram"))

# combine and write
ss <-
  ss %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(donor_id = donor_id,
                mutations = paste0(wd, "/out/nf-resolveome/PD63118/mutations.tsv"))
ss %>%
  readr::write_tsv(paste0(out_dir, "/samplesheet.tsv"))
