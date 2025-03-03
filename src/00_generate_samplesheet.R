#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
data_dir <- paste0(wd, "/data/resolveome/")

# read in manifests
manifest <-
  list.files("data/manifest/", pattern = ".xlsx$", full.names = TRUE) %>%
  {setNames(., basename(.) %>% tools::file_path_sans_ext())} %>%
  purrr::map(function(file) {
    file %>%
      readxl::read_excel(skip = 8) %>%
      janitor::clean_names()
  })

# run49686_lane3 (7761stdy_manifest_25650_031024)
# filter samples and fix supplier sample names
manifest[["7761stdy_manifest_25650_031024"]] <-
  manifest[["7761stdy_manifest_25650_031024"]] %>%
  dplyr::filter(grepl("^Hashimoto", supplier_sample_name)) %>%
  dplyr::mutate(
    supplier_sample_name = paste0("PD63118_HYB_",
                                  gsub(".* - ", "", supplier_sample_name)))

# run49686_lane4-5 (7761stdy_manifest_25654_041024)
# fix supplier sample names
manifest[["7761stdy_manifest_25654_041024"]] <-
  manifest[["7761stdy_manifest_25654_041024"]] %>%
  dplyr::mutate(
    supplier_sample_name = paste0("PD63118_", gsub(" - ", "_",
                                                   supplier_sample_name)))

# run 49901_lane8 (7894stdy_manifest_26015_211124_RNA)
# remove control (!= "1 cell") and low-yield (quant < 2) wells
# also, there is a typo in the manifest, where 7894STDY15290419 and
# 7894STDY15290420 are both assigned as PD63118_RNA_F10, while no sample is
# assigned to PD63118_RNA_G10. i will fix this by assigning 7894STDY15290420 to
# PD63118_RNA_G10.
pre_pcr_quants <-
  "data/plate_layout/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout_dna_pre_pcr_quants.tsv" %>%
  readr::read_tsv() %>%
  tidyr::pivot_longer(cols = -`...1`, values_to = "well_quant") %>%
  dplyr::mutate(well = paste0(`...1`, name))
plate_layout <-
  "data/plate_layout/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout_plate_layout.tsv" %>%
  readr::read_tsv() %>%
  tidyr::pivot_longer(cols = -`...1`, values_to = "well_content") %>%
  dplyr::mutate(well = paste0(`...1`, name))
pass_samples <-
  dplyr::full_join(pre_pcr_quants, plate_layout) %>%
  dplyr::filter(well_quant >= 2, well_content == "1 cell") %>%
  dplyr::mutate(supplier_sample_name = paste0("PD63118_RNA_", well)) %>%
  dplyr::pull(supplier_sample_name)
manifest[["7894stdy_manifest_26015_211124_RNA"]] <-
  manifest[["7894stdy_manifest_26015_211124_RNA"]] %>%
  dplyr::mutate(
    supplier_sample_name = ifelse(sanger_sample_id == "7894STDY15290420",
                                  "PD63118_RNA_G10", supplier_sample_name)) %>%
  dplyr::filter(supplier_sample_name %in% pass_samples)

# combine standardised manifests
manifest <-
  manifest %>%
  dplyr::bind_rows(.id = "manifest_file")

# read in sequencescape pool files
seqscape <-
  list.files("data/sequencescape/", pattern = "tsv$", full.names = TRUE) %>%
  {purrr::set_names(., basename(.) %>% tools::file_path_sans_ext())} %>%
  purrr::map(readr::read_tsv) %>%
  dplyr::bind_rows(.id = "id") %>%
  janitor::clean_names() %>%
  tidyr::separate_wider_delim("id", delim = "_", names = c("run_id", "lane"))

# hard-code the run_id-to-plate conversion
run_id_to_plate <-
  c("49686" = 1, "49882" = 3, "49900" = 3, "49901" = 3, "50072" = 3)

# combine manifest and sequencescape
# generate cell_id - plate{plate}_well{well}
# /seq/illumina/runs/49/49686/lane4-5/plex1/49686_4-5#1.cram
ss <-
  dplyr::inner_join(manifest, seqscape) %>%
  tidyr::separate_wider_delim("supplier_sample_name", delim = "_",
                              names = c("pdid", "seq_type", "well"),
                              cols_remove = FALSE) %>%
  dplyr::mutate(
    # specify all HYB sequencing is from DNA
    seq_type = ifelse(seq_type == "HYB", "DNAHYB", seq_type) %>% tolower(),
    lane_tmp = ifelse(lane == "lane1-8", "", lane),
    lane_n_tmp = ifelse(lane == "lane1-8", "", paste0("_", gsub("lane", "",
                                                                lane))),
    bam = paste0("/seq/illumina/runs/", substr(run_id, 1, 2), "/", run_id, "/",
                 lane_tmp, "/plex", npg_aliquot_index, "/", run_id, lane_n_tmp,
                 "#", npg_aliquot_index, ".cram"),
    plate = run_id_to_plate[run_id]) %>%
  dplyr::transmute(
    cell_id = paste0("plate", plate, "_well", well),
    id = paste0(cell_id, "_", seq_type, "_run", run_id),
    plate, well, seq_type, run_id, lane, plex_n = npg_aliquot_index,
    study_id = gsub("STDY.*", "", sanger_sample_id),
    donor_id = donor_id_required_for_ega,
    sanger_sample_id, supplier_sample_name, manifest_file,
    bam)

# write samplesheets
ss %>%
  readr::write_csv(paste0(data_dir, "/samplesheet_irods.csv"))
ss_local <-
  ss %>%
  dplyr::mutate(bam = paste0(data_dir, "/", donor_id, "/", id, "/bam/", id,
                             ".bam"))
ss_local %>%
  readr::write_csv(paste0(data_dir, "/samplesheet_local.csv"))

# write a samplesheet updating the merged bams
ss_merged <-
  ss_local %>%
  dplyr::filter(run_id %in% c(49901, 50072), seq_type == "rna") %>%
  dplyr::mutate(
    id = paste0(cell_id, "_", seq_type, "_merged"),
    bam = paste0(wd, "/out/merge_bams/", donor_id, "/", id, "/", "/bam/", id,
                 ".bam")) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(dplyr::across(everything(), ~ paste(unique(.),
                                                       collapse = ",")))
ss_unmerged <-
  ss_local %>%
  dplyr::filter(!(run_id %in% c(49901, 50072) & seq_type == "rna")) %>%
  dplyr::mutate(dplyr::across(everything(), as.character))
ss_w_merged <- dplyr::bind_rows(ss_merged, ss_unmerged)
ss_w_merged %>%
  readr::write_csv(paste0(data_dir, "/samplesheet_local_merged.csv"))

# write samplesheet for testing (3 samples from each run)
ss_local %>%
  dplyr::group_by(seq_type, run_id, lane) %>%
  dplyr::filter(dplyr::row_number() < 4) %>%
  readr::write_csv("out/test/samplesheet.csv")
