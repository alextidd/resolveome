# libraries
library(magrittr)

# dirs
donor_id <- "PD63118"
out_dir <- paste0("out/nf-resolveome/", donor_id, "/")
dir.create(out_dir)
wd <- getwd()

# function: define mutation type based on ref and alt, split up mnvs and dnvs
type_mutations <- function(df) {
  typed_df <-
    df %>%
    dplyr::mutate(type = dplyr::case_when(
      nchar(ref) == 1 & nchar(mut) == 1 ~ "snv",
      nchar(ref) == 1 & nchar(mut) > 1 ~ "ins",
      nchar(ref) > 1 & nchar(mut) == 1 ~ "del",
      nchar(ref) == 2 & nchar(mut) == 2 ~ "dnv",
      nchar(ref) > 1 & nchar(mut) > 1 ~ "mnv"
    ))

  # expand dnv/mnv mutations to all positions
  if (any(typed_df$type %in% c("dnv", "mnv"))) {
    typed_df_dnv_mnv <-
      typed_df %>%
      dplyr::filter(type %in% c("dnv", "mnv")) %>%
      dplyr::mutate(mut_id = paste(chr, pos, ref, mut, type)) %>%
      dplyr::group_by(dplyr::across(-c(pos, mut, ref))) %>%
      dplyr::reframe(pos = pos:(pos + nchar(mut) - 1),
                    ref = strsplit(ref, split = "") %>% unlist(),
                    mut = strsplit(mut, split = "") %>% unlist()) %>%
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
caveman_snps <-
  "/nfs/cancer_ref01/nst_links/live/3438/PD63118b_lo0009/PD63118b_lo0009.caveman_c.snps.vcf.gz" %>%
  readr::read_tsv(comment = "##") %>%
  dplyr::mutate(
    DP = strsplit(INFO, ";") %>% purrr::map_chr(~ .x[grepl("^DP=", .x)]) %>%
      strsplit("=") %>% purrr::map_chr(~ .x[2]) %>% as.integer(),
    VAF = gsub(".*:", "", TUMOUR) %>% as.numeric()) %>%
  dplyr::filter(DP > 50, 0.3 < VAF, VAF < 0.7) %>%
  # get those at common snp sites
  dplyr::inner_join(common_snps) %>%
  dplyr::transmute(donor_id = "PD63118", chr = `#CHROM`, pos = POS, ref = REF,
                   mut = ALT) %>%
  # type the mutations
  type_mutations() %>%
  dplyr::distinct()

# nanoseq samplesheet
ss_muts <- readr::read_tsv("data/nanoseq/samplesheet.tsv")

# get mutations from exome and targeted nanoseq
nanoseq_muts <-
  ss_muts %>%
  purrr::pmap(function(muts_file, seq_type) {
    readr::read_tsv(muts_file) %>%
      dplyr::mutate(seq_type = seq_type)
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(donor_id = stringr::str_sub(sampleID, 1, 7)) %>%
  dplyr::filter(donor_id == curr_donor_id) %>%
  dplyr::distinct(donor_id, chr, pos, ref, mut) %>%
  # type the mutations
  type_mutations()

# write mutations
muts_file <- file.path(wd, out_dir, donor_id, "mutations.tsv")
muts_and_snps <-
  list("caveman_snps" = caveman_snps,
       "nanoseq_mutations" = nanoseq_muts) %>%
  dplyr::bind_rows(.id = "source")
muts_and_snps %>%
  readr::write_tsv(muts_file)

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

# combine and write
ss <-
  ss %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(donor_id = paste0(donor_id, "_", run),
                id = paste0("plex", plex_n),
                mutations = muts_file)
ss %>% readr::write_tsv(paste0(out_dir, "/samplesheet.tsv"))
