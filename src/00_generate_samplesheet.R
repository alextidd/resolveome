#!/usr/bin/env Rscript

# libraries
library(magrittr)

# dirs
wd <- getwd()
data_dir <- paste0(wd, "/data/resolveome/")

# read in manifests
manifests <-
  list.files("data/manifest/", pattern = ".xlsx$", full.names = TRUE) %>%
  {setNames(., basename(.) %>% tools::file_path_sans_ext())} %>%
  purrr::map(function(file) {
    file %>%
      readxl::read_excel(skip = 8) %>%
      janitor::clean_names()
  })

# run49686_lane3 (7761stdy_manifest_25650_031024)
# filter samples and fix supplier sample names
manifests[["7761stdy_manifest_25650_031024"]] <-
  manifests[["7761stdy_manifest_25650_031024"]] %>%
  dplyr::filter(grepl("^Hashimoto", supplier_sample_name)) %>%
  dplyr::mutate(
    supplier_sample_name = paste0("PD63118_HYB_",
                                  gsub(".* - ", "", supplier_sample_name)))

# run49686_lane4-5 (7761stdy_manifest_25654_041024)
# fix supplier sample names
manifests[["7761stdy_manifest_25654_041024"]] <-
  manifests[["7761stdy_manifest_25654_041024"]] %>%
  dplyr::mutate(
    supplier_sample_name = paste0("PD63118_", gsub(" - ", "_",
                                                   supplier_sample_name)))

# run50227_lane7 (7761stdy_manifest_26617_Hyb_manifest_run50227_lane7_corrected)
# fix supplier sample names
manifests[["7761stdy_manifest_26617_Hyb_manifest_run50227_lane7_corrected"]] <-
  manifests[["7761stdy_manifest_26617_Hyb_manifest_run50227_lane7_corrected"]] %>%
  dplyr::mutate(
    supplier_sample_name = gsub("_P10", "", supplier_sample_name))

# run50227_lane8 (7761stdy_manifest_26618_RNA_manifest_run50227_lane8_corrected)
# fix supplier sample names
manifests[["7761stdy_manifest_26618_RNA_manifest_run50227_lane8_corrected"]] <-
  manifests[["7761stdy_manifest_26618_RNA_manifest_run50227_lane8_corrected"]] %>%
  dplyr::mutate(
    supplier_sample_name = gsub("_P10", "", supplier_sample_name))

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
manifests[["7894stdy_manifest_26015_211124_RNA"]] <-
  manifests[["7894stdy_manifest_26015_211124_RNA"]] %>%
  dplyr::mutate(
    supplier_sample_name = ifelse(sanger_sample_id == "7894STDY15290420",
                                  "PD63118_RNA_G10", supplier_sample_name)) %>%
  dplyr::filter(supplier_sample_name %in% pass_samples)

# combine standardised manifests
manifest <-
  manifests %>%
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
  c("49686" = 1,
    "49882" = 3, "49900" = 3, "49901" = 3, "50072" = 3,
    "50227" = 10)

# ss irods - combine manifest and sequencescape
# generate cell_id - plate{plate}_well{well}
# generate id      - plate{plate}_well{well}_{seq_type}_run{run_id}
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
    # merge the two duplicate RNA runs
    id = dplyr::case_when(
      run_id %in% c(49901, 50072) & seq_type == "rna" ~
        paste0(cell_id, "_", seq_type, "_merged"),
      TRUE ~ paste0(cell_id, "_", seq_type, "_run", run_id)),
    plate, well, seq_type, run_id, lane, plex_n = npg_aliquot_index,
    study_id = gsub("STDY.*", "", sanger_sample_id),
    donor_id = donor_id_required_for_ega,
    sanger_sample_id, supplier_sample_name, manifest_file,
    bam)

# write ss irods
ss %>%
  readr::write_csv(paste0(data_dir, "/samplesheet_irods.csv"))

# ss local (with merged RNA BAMs collapsed)
ss_local <-
  ss %>%
  dplyr::mutate(bam = paste0(data_dir, "/", donor_id, "/", id, "/bam/", id,
                             ".bam")) %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(dplyr::across(everything(), ~ paste(unique(.),
                                                       collapse = "|")))

# write ss local
ss_local %>%
  readr::write_csv(paste0(data_dir, "/samplesheet_local.csv"))

# write samplesheet for testing (3 samples from each run)
ss_local %>%
  dplyr::group_by(seq_type, run_id, lane) %>%
  dplyr::filter(dplyr::row_number() < 4) %>%
  readr::write_csv("out/test/samplesheet.csv")

# wrangle plate10 manual inspection
plate10_man_insp <-
  readr::read_tsv("data/manual_inspection/2025-03-28_PTA_PD63118_Plate10.tsv")

# add pre pcr quants
dnahyb_pre_pcr_quants <-
  "data/plate_layout/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout_dna_pre_pcr_quants.tsv" %>%
  readr::read_tsv() %>%
  tidyr::pivot_longer(cols = -`...1`, values_to = "DNA_PrePCR_conc") %>%
  dplyr::transmute(well = paste0(`...1`, name), DNA_PrePCR_conc)
rna_pre_pcr_quants <-
  "data/plate_layout/2024-11-13_Hashimoto_PD63118_Plate3_PlateLayout_rna_pre_pcr_quants.tsv" %>%
  readr::read_tsv() %>%
  tidyr::pivot_longer(cols = -`...1`, values_to = "RNA_PrePCR_conc") %>%
  dplyr::transmute(well = paste0(`...1`, name), RNA_PrePCR_conc)

# combine with ids
plate10_man_insp %>%
  dplyr::transmute(run_id = as.character(run_id), plex = paste0("plex", plex),
                   n_cells, `1p_loh`, TNFRSF14_mut) %>%
  dplyr::inner_join(
    ss %>%
      dplyr::filter(seq_type == "dnahyb") %>%
      dplyr::mutate(plex = paste0("plex", plex_n))) %>% 
  dplyr::left_join(dnahyb_pre_pcr_quants) %>%
  dplyr::left_join(rna_pre_pcr_quants) %>%
  dplyr::transmute(
    id, run_id, n_cells, plex, cell_id, seq_type, DNA_PrePCR_conc, RNA_PrePCR_conc, suspected_doublet = "",
    doublet_rationale = "", celltype_SHM = "", celltype_VDJ_recomb = "", class_switch_recombination_CSR = "", chr_dropout = "",
    `1p_loh`, TNFRSF14_mut, TNFRSF14_mut_VAF = "", TNFRSF14_mut_in_NanoSeq_data = "", productive_heavy_chain = "",
    heavy_chain_CDR3_nt_IgBLAST = "", heavy_chain_CDR3_aa_IgBLAST = "", unproductive_heavy_chain = "", heavy_chain_isotype = "",
    productive_light_chain = "", identity_of_BCR_light_chain = "", productive_TRB_chain = "", productive_TRA_chain = "",
    TRG_rearrangement = "", notes = "") %>%
  readr::write_tsv("data/manual_inspection/2025-03-28_PTA_PD63118_Plate10_reformatted.tsv")
