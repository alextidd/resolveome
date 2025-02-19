#!/usr/bin/env Rscript

# libraries
library(magrittr)

# samplesheet
ss_bams <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type %in% c("dna", "dnahyb"))
ss_muts <- readr::read_tsv("data/nanoseq/samplesheet.tsv")

# dirs
donor_id_i <- "PD63118"
out_dir <- "out/nf-resolveome/"
dir.create(out_dir, showWarnings = FALSE)
wd <- getwd()

# function: define mutation type based on ref and alt, split up mnvs and dnvs
type_mutations <- function(df) {
  typed_df <-
    df %>%
    dplyr::mutate(type = dplyr::case_when(
      nchar(ref) == 1 & nchar(alt) == 1 ~ "snv",
      nchar(ref) == 1 & nchar(alt) > 1 ~ "ins",
      nchar(ref) > 1 & nchar(alt) == 1 ~ "del",
      nchar(ref) == 2 & nchar(alt) == 2 ~ "dnv",
      nchar(ref) > 1 & nchar(alt) > 1 ~ "mnv"
    ))

  # expand dnv/mnv mutations to all positions
  if (any(typed_df$type %in% c("dnv", "mnv"))) {
    typed_df_dnv_mnv <-
      typed_df %>%
      dplyr::filter(type %in% c("dnv", "mnv")) %>%
      dplyr::mutate(mut_id = paste(chr, pos, ref, alt, type)) %>%
      dplyr::group_by(dplyr::across(-c(pos, alt, ref))) %>%
      dplyr::reframe(pos = pos:(pos + nchar(alt) - 1),
                     ref = strsplit(ref, split = "") %>% unlist(),
                     alt = strsplit(alt, split = "") %>% unlist()) %>%
      dplyr::select(-mut_id)
    typed_df <-
      typed_df %>%
      dplyr::filter(!(type %in% c("dnv", "mnv"))) %>%
      dplyr::bind_rows(typed_df_dnv_mnv)
  }

  # return
  typed_df
}

# get common snp sites
common_snps <-
  gzfile("/lustre/scratch125/casm/team268im/fa8/117/NOVASEQX_MASKS/NEW_MASKS_wNSX_OCT2024/GRCh37_WGNS/SNP_GRCh37.wgns.bed.gz") %>%
  readr::read_tsv(col_names = c("#CHROM", "START", "POS")) %>%
  dplyr::select(`#CHROM`, POS)

# get a caveman snp file from a sample with high coverage
# extract common snps that are heterozygous
# 0.3 < VAF < 0.7 and DP > 50
# PD63118b_lo0044 has the highest coverage at 68X according to picard
caveman_snps <-
  "/nfs/cancer_ref01/nst_links/live/3464/PD63118b_lo0044/PD63118b_lo0044.caveman_c.snps.vcf.gz" %>%
  readr::read_tsv(comment = "##") %>%
  dplyr::mutate(
    DP = strsplit(INFO, ";") %>% purrr::map_chr(~ .x[grepl("^DP=", .x)]) %>%
      strsplit("=") %>% purrr::map_chr(~ .x[2]) %>% as.integer(),
    VAF = gsub(".*:", "", TUMOUR) %>% as.numeric()) %>%
  dplyr::filter(DP > 50, 0.3 < VAF, VAF < 0.7) %>%
  # get those at common snp sites
  dplyr::inner_join(common_snps) %>%
  dplyr::transmute(donor_id = donor_id_i, chr = `#CHROM`, pos = POS, ref = REF,
                   alt = ALT) %>%
  # type the mutations
  type_mutations() %>%
  dplyr::distinct()

# get mutations from exome and targeted nanoseq
nanoseq_muts <-
  ss_muts %>%
  purrr::pmap(function(muts_file, seq_type) {
    readr::read_tsv(muts_file) %>%
      dplyr::mutate(seq_type = seq_type)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(donor_id = stringr::str_sub(sampleID, 1, 7)) %>%
  dplyr::filter(donor_id == donor_id_i) %>%
  dplyr::rename(alt = mut) %>%
  dplyr::distinct(donor_id, chr, pos, ref, alt) %>%
  # type the mutations
  type_mutations()

# write mutations
muts_file <- file.path(wd, out_dir, donor_id_i, "mutations.tsv")
muts_and_snps <-
  list("caveman_snps" = caveman_snps,
       "nanoseq_mutations" = nanoseq_muts) %>%
  dplyr::bind_rows(.id = "source")
muts_and_snps %>%
  readr::write_tsv(muts_file)

# combine and write
ss <-
  ss_bams %>% dplyr::mutate(mutations = muts_file)
ss %>% readr::write_csv(file.path(out_dir, "samplesheet.csv"))
