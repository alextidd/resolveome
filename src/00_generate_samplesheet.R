#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
data_dir <- "data/resolveome/DNA/"
dir.create(data_dir)

# initiate samplesheet list
ss <- list()

# curr pdid
donor_id_i <- "PD63118"

# get plex conversion between DNA and RNA for plate3 (49900/499001)
file_base_name <- "2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout"
dfs <-
  paste0("data/manifest/", file_base_name,
         c("_dna_pre_pcr_quants.tsv", "_rna_pre_pcr_quants.tsv",
           "_plate_layout.tsv")) %>%
  {setNames(., basename(.) %>% gsub(paste0(file_base_name, "_"), "", .) %>% tools::file_path_sans_ext())} %>%
  purrr::map(function(df) {
    readr::read_tsv(df) %>%
      dplyr::rename(row = `...1`) %>%
      tidyr::pivot_longer(cols = as.character(1:12),
                          names_to = "column") %>%
      dplyr::mutate(column = as.numeric(column),
                    well_id = paste0(row, column)) %>%
      dplyr::arrange(column, row)
  })
plate3 <-
  dfs$plate_layout %>%
  dplyr::rename(contents = value) %>%
  dplyr::left_join(dfs$dna_pre_pcr_quants %>% dplyr::rename(dna_quant = value)) %>%
  dplyr::left_join(dfs$rna_pre_pcr_quants %>% dplyr::rename(rna_quant = value)) %>%
  dplyr::mutate(
    # DNA plex goes up by rows and then columns, with empty wells and low quant
    # wells (<2) excluded
    dna_plex_n = dplyr::case_when(contents == "1 cell" & dna_quant > 2 ~ 1,
                                  TRUE ~ NA),
    dna_plex_n = cumsum(dplyr::coalesce(dna_plex_n, 0)) + dna_plex_n * 0,
    # RNA plex goes up by rows and then columns, with empty wells included
    rna_plex_n = dplyr::case_when(contents %in% c("EMPTY", "1 cell") ~ 1,
                                  TRUE ~ NA),
    rna_plex_n = cumsum(dplyr::coalesce(rna_plex_n, 0)) + rna_plex_n * 0,
    dna_run = "49900", rna_run = "49901")

# save
plate3 %>%
  readr::write_tsv("data/manifest/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout.tsv")

# save irods samplesheet
plate3 %>%
  dplyr::filter(!is.na(rna_plex_n), !is.na(dna_plex_n)) %>%
  dplyr::transmute(
    id = paste0("plex", dna_plex_n), lane = 8,
    run = rna_run, donor_id = paste0("PD63118_", run),
    bam = paste0("/seq/illumina/runs/", substr(run, 1, 2), "/", run, "/lane",
                 lane, "/plex", rna_plex_n, "/", run, "_", lane, "#",
                 rna_plex_n, ".cram")) %>%
  readr::write_csv("data/resolveome/RNA/samplesheet_irods.csv")

# samplesheet 49686 (19 cells)
ss[["49686"]] <-
  tibble::tibble(
    run = "49686",
    lane_n = "4-5",
    plex_n = seq(1, 19)) %>%
  dplyr::transmute(
    run, plex_n, lane_n,
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
    bam = paste0("/seq/illumina/runs/", substr(run, 1, 2), "/", run, "/lane",
                 lane_n, "/plex", plex_n, "/", run, "_", lane_n, "#", plex_n,
                 ".cram"))

# combine
ss <-
  ss %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(donor_id = paste0(donor_id_i, "_", run),
                id = paste0("plex", plex_n))

# write irods samplesheet
ss %>% readr::write_csv(file.path(data_dir, "/samplesheet_irods.csv"))

# write local samplesheet
ss %>%
  dplyr::mutate(bam = file.path(wd, data_dir, donor_id, id, "bam",
                                paste0(id, ".bam"))) %>%
  readr::write_csv(file.path(data_dir, "/samplesheet_local.csv"))
