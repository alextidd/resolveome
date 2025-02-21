# libraries
library(magrittr)
library(ImmuneLENS)

# point to samtools module
samtools_dir <- "/software/spack_environments/default/00/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-13.1.0/samtools-1.19.2-2b4gy7kosouxuikfcv645lg23ygnk2qu/bin/"
Sys.setenv(PATH = paste(samtools_dir, Sys.getenv("PATH"), sep = ":"))

# dirs
out_dir <- "out/immunelens/"
dir.create(out_dir, showWarnings = FALSE)

# samplesheet
ss <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type == "dna")

# test bams
# B cell: plate3_wellF2_dna_run49882
# T cell: plate3_wellD2_dna_run49882
test_bams <-
  ss %>%
  dplyr::filter(id %in% c("plate3_wellF2_dna_run49882",
                          "plate3_wellD2_dna_run49882")) %>%
  dplyr::pull(bam) %>%
  {purrr::set_names(., tools::file_path_sans_ext(basename(.)))}

pdf("test.pdf")
results <-
  purrr::map(test_bams, function(test_bam) {

  # test_id <- "plate3_wellA2_dna_run49882" ; test_bam <- paste0("/lustre/scratch125/casm/team268im/at31/resolveome/data/resolveome/PD63118/", test_id, "/bam/", test_id, ".bam")
  # get coverage
  TCRA.cov <- getCovFromBam_WGS(bamPath = test_bam, outPath = out_dir,
                                vdj.gene = "TCRA", hg19_or_38 = "hg19")
  TCRB.cov <- getCovFromBam_WGS(bamPath = test_bam, outPath = out_dir,
                                vdj.gene = "TCRB", hg19_or_38 = "hg19")
  TCRG.cov <- getCovFromBam_WGS(bamPath = test_bam, outPath = out_dir,
                                vdj.gene = "TCRG", hg19_or_38 = "hg19")
  IGH.cov <- getCovFromBam_WGS(bamPath = test_bam, outPath = out_dir,
                               vdj.gene = "IGH", hg19_or_38 = "hg19")

  # load coverage
  TCRA.df <- loadCov(TCRA.cov)
  TCRB.df <- loadCov(TCRB.cov)
  TCRG.df <- loadCov(TCRG.cov)

  # run T cell ExTRECT for TCRA, TCRB, TCRG using new WGS version
  TCRA.out <- runImmuneLENS(vdj.region.df = TCRA.df, vdj.gene = "TCRA",
                            hg19_or_38 = "hg19")
  TCRB.out <- runImmuneLENS(vdj.region.df = TCRB.df, vdj.gene = "TCRB",
                            hg19_or_38 = "hg19", median.thresh = 0)
  TCRG.out <- runImmuneLENS(vdj.region.df = TCRG.df, vdj.gene = "TCRG",
                            hg19_or_38 = "hg19")

  # plot
  plotImmuneLENS(vdj.region.df = TCRA.df, vdj.gene = "TCRA",
                hg19_or_38 = "hg38", ylims = c(-1.5, 1)) %>% print()
  plotImmuneLENS(vdj.region.df = TCRB.df, vdj.gene = "TCRB",
                hg19_or_38 = "hg38", ylims = c(-1.5, 1)) %>% print()
  plotImmuneLENS(vdj.region.df = TCRG.df, vdj.gene = "TCRG",
                hg19_or_38 = "hg38", ylims = c(-1.5, 1)) %>% print()

  # return
  list(TCRA = TCRA.out, TCRB = TCRB.out, TCRG = TCRG.out)

})
dev.off()
