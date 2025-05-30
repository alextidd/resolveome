# cd /lustre/scratch125/casm/team268im/at31/resolveome ; bsub -q basement -M32000 -R 'span[hosts=1] select[mem>32000] rusage[mem=32000]' -J 06_ImmuneLENS_run -o log/%J_06_ImmuneLENS_run.out -e log/%J_06_ImmuneLENS_run.err 'module load ISG/rocker/rver/4.4.0; export R_LIBS_USER=$HOME/R-tmp-4.4; Rscript src/06_ImmuneLENS_run.R'

# libraries
library(magrittr)
library(ImmuneLENS)

# point to samtools module
samtools_dir <- "/software/spack_environments/default/00/opt/spack/linux-ubuntu22.04-x86_64_v3/gcc-13.1.0/samtools-1.19.2-2b4gy7kosouxuikfcv645lg23ygnk2qu/bin/"
Sys.setenv(PATH = paste(samtools_dir, Sys.getenv("PATH"), sep = ":"))

# dirs
out_dir <- "out/immunelens/"
cov_dir <- file.path(out_dir, "cov/")
cov_plots_dir <- file.path(out_dir, "cov_plots/")
extrect_dir <- file.path(out_dir, "extrect/")
dir.create(out_dir, showWarnings = FALSE)
dir.create(cov_dir, showWarnings = FALSE)
dir.create(cov_plots_dir, showWarnings = FALSE)
dir.create(extrect_dir, showWarnings = FALSE)

# samplesheet
ss <-
  readr::read_csv("data/resolveome/samplesheet_local.csv") %>%
  dplyr::filter(seq_type == "dna") %>%
  tail(56)

# get genes
genes <- c("TCRA", "TCRB", "TCRG", "IGH")

# get scores
scores <-
  purrr::map2(ss$id, ss$bam, function(id_i, bam) {

    print(paste(id_i, "-", bam))
    message(paste(id_i, "-", bam))

    genes %>%
      purrr::map(function(gene) {

        # id_i <- ss$id[1] ; bam <- ss$bam[1] ; gene <- genes[1]
        # id_i <- "plate1_wellC2_dna_run49686" ; gene <- "TCRG" ; bam <- dplyr::filter(ss, id == id_i)$bam
        # id_i <- "plate3_wellC11_dna_run49882" ; gene <- "TCRG" ; bam <- "/lustre/scratch125/casm/team268im/at31/resolveome/data/resolveome//PD63118/plate3_wellC11_dna_run49882/bam/plate3_wellC11_dna_run49882.bam"

        print(gene)
        message(gene)

        # set celltype
        ct <- ifelse(gene == "IGH", "bcell", "tcell")

        print("get coverage")
        cov_file <-
          getCovFromBam_WGS(bamPath = bam, outPath = cov_dir, vdj.gene = gene,
                            hg19_or_38 = "hg19")

        print("load coverage")
        cov <- suppressMessages(loadCov(cov_file))

        print("plot coverage")
        p <-
          tryCatch(
            plotImmuneLENS(vdj.region.df = cov, vdj.gene = gene,
                           median.thresh = 0, hg19_or_38 = "hg19",
                           sample_name = id_i),
            error = function(e) {
              message("Plotting error caught: ", conditionMessage(e))
              NULL
            }
          )
        if (!is.null(p)) {
          png(paste0(cov_plots_dir, id_i, "_", gene, "_cov.png"),
              width = 10, height = 5, units = "in", res = 300)
          print(p)
          dev.off()
        }

        print("run T cell ExTRECT, catch errors")
        extrect <-
          tryCatch(
            runImmuneLENS(vdj.region.df = cov, vdj.gene = gene,
                          hg19_or_38 = "hg19", median.thresh = 0),
            error = function(e) {
              message("Error caught: ", conditionMessage(e))
              list(tibble::tibble(name = paste0(gene, ".", ct, ".fraction"),
                                  value = NA) %>% tidyr::pivot_wider())
            })

        print("save extrect results")
        extrect %>%
          saveRDS(file = file.path(extrect_dir, paste0(id_i, ".", gene, ".rds")))

        print("DONE!")
        tibble::tibble(
          cell_fraction = extrect[[1]][, paste0(gene, ".", ct, ".fraction"),
                                       drop = TRUE],
          celltype = ct, gene = gene, id = id_i)

    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()

# save scores
scores %>% readr::write_tsv(file.path(out_dir, "scores.tsv"))