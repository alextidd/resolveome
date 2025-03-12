library(dplyr)
library(magrittr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

cutoff_coverage <- args[1] %>% as.numeric()
cutoff_variant <- args[2] %>% as.numeric()
cutofof_prev_na <- args[3] %>% as.numeric()
cutoff_prev_var <- args[4] %>% as.character()
mchrs <- args[5] %>% as.character()

mfiles <- list.files(pattern = "res_grouplevel*") %>% head(30)

msamples <- mfiles %>%
  gsub(pattern = ".*_sample_", replacement = "") %>%
  gsub(pattern = "_file.*", replacement = "") %>%
  unique() %>%
  sort()

moutfile_nr <- "Tab_NR_somatic_snps.tsv"
moutfile_nv <- "Tab_NV_somatic_snps.tsv"

for (mchr in mchrs) {
  sub_files <- mfiles %>% grep(pattern = paste0("_chr_", mchr, "_sample"), value = TRUE)

  df_pileup <- NULL

  i <- 0

  total_files <- length(sub_files)


  for (mfile in sub_files) {
    i <- i + 1
    cat("Working on chr: ", mchr, "\t file: ", mfile, "\t ", i, " out of", total_files, "\n")

    mname <- mfile %>%
      gsub(pattern = ".*_sample_", replacement = "") %>%
      gsub(pattern = "_file.*", replacement = "")

    df_temp <- data.table::fread(mfile, sep = "\t", quote = "") %>%
      as.data.frame() %>%
      dplyr::mutate(.data = ., UId = paste0(CHROM, "_", POS, "_", REF, "_", ALT)) %>%
      dplyr::mutate(.data = ., UIdPos = paste0(CHROM, "_", POS))


    df_pileup <- with(df_temp, order(POS)) %>%
      df_temp[., ] %>%
      dplyr::select(.data = ., c("UId", "UIdPos", "DEPTH", "NSUPVAR", "NSUPREF")) %>%
      dplyr::mutate(.data = ., SampleId = mname) %>%
      dplyr::mutate(.data = ., UIdPosSampleId = paste0(UIdPos, "_", SampleId)) %>%
      rbind(df_pileup, .)

    rm(df_temp)
  }

  df_depth <- df_pileup[, c("UIdPos", "DEPTH", "SampleId", "UIdPosSampleId")] %>% unique()



  df_nv <- df_pileup$UId %>%
    grep(pattern = "REF", invert = TRUE) %>%
    df_pileup[., ] %>%
    subset(NSUPVAR != 0)

  ids_prob <- df_depth[, c("UIdPosSampleId")] %>%
    table() %>%
    data.frame() %>%
    subset(Freq > 1) %$% .

  # Recalculate depth for this ones
  ids_good <- df_depth[, c("UIdPosSampleId")] %>%
    table() %>%
    data.frame() %>%
    subset(Freq == 1) %$% .


  # Subset properly contained
  df_depth_clean <- df_depth %>%
    dplyr::filter(.data = ., UIdPosSampleId %in% ids_good) %>%
    droplevels()

  # Solve missing ones using pileup information
  df_prob <- df_pileup %>%
    dplyr::filter(.data = ., UIdPosSampleId %in% ids_prob) %>%
    droplevels()

  df_a <- aggregate(NSUPVAR ~ UIdPos + SampleId, df_prob, sum)
  df_b <- aggregate(NSUPREF ~ UIdPos + SampleId, df_prob, max)

  merged <- merge(df_a, df_b, by = c("UIdPos", "SampleId"))

  merged$DEPTH <- merged$NSUPVAR + merged$NSUPREF

  df_depth <- rbind(df_depth_clean[, c("UIdPos", "SampleId", "DEPTH")], merged[, c("UIdPos", "SampleId", "DEPTH")])

  # Create TABS
  Tab_NR <- reshape2::acast(df_depth, UIdPos ~ SampleId, value.var = "DEPTH", fill = 0)

  Tab_NV <- reshape2::acast(df_nv, UId ~ SampleId, value.var = "NSUPVAR", fill = 0)

  # Now create the missing part
  mpos <- rownames(Tab_NV) %>%
    strsplit(split = "_") %>%
    lapply(function(x) paste0(x[1], "_", x[2])) %>%
    unlist()

  Tab_NR <- match(mpos, rownames(Tab_NR)) %>% Tab_NR[., ]

  rownames(Tab_NR) <- rownames(Tab_NV)

  Tab_NR <- match(msamples, colnames(Tab_NR)) %>% Tab_NR[, .]

  Tab_NV <- match(msamples, colnames(Tab_NV)) %>% Tab_NV[, .]

  # Remove NA_NA from appending
  Tab_NR <- rownames(Tab_NR) %>%
    grep(pattern = "NA_NA", invert = TRUE) %>%
    Tab_NR[., ]

  Tab_NV <- rownames(Tab_NV) %>%
    grep(pattern = "NA_NA", invert = TRUE) %>%
    Tab_NV[., ]

  if (mchr == "chr1") {
    # Write matrices
    mheader <- c("", colnames(Tab_NR))
    write.table(x = mheader %>% t(), file = moutfile_nr, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
    mheader <- c("", colnames(Tab_NV))
    write.table(x = mheader %>% t(), file = moutfile_nv, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  }
  write.table(x = Tab_NR, file = moutfile_nr, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
  write.table(x = Tab_NV, file = moutfile_nv, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
}

### Second part to pre-filter
moutfile_nr <- "Tab_NR_somatic_snps_filtered.tsv"
moutfile_nv <- "Tab_NV_somatic_snps_filtered.tsv"

if (mchrs[1] == "chr1") {
  Tab_nr <- read.table("Tab_NR_somatic_snps.tsv", sep = "\t", header = TRUE, row.names = 1, comment.char = "", check.names = FALSE) %>%
    as.matrix()

  Tab_nv <- read.table("Tab_NV_somatic_snps.tsv", sep = "\t", header = TRUE, row.names = 1, comment.char = "", check.names = FALSE) %>%
    as.matrix()
} else {
  Tab_nr <- read.table("Tab_NR_somatic_snps.tsv", sep = "\t", header = FALSE, row.names = 1, comment.char = "", check.names = FALSE) %>%
    as.matrix()

  Tab_nv <- read.table("Tab_NV_somatic_snps.tsv", sep = "\t", header = FALSE, row.names = 1, comment.char = "", check.names = FALSE) %>%
    as.matrix()
}



Tab_nr_copy <- Tab_nr
Tab_nv_copy <- Tab_nv

Tab_nr_copy[Tab_nr_copy < cutoff_coverage] <- NA

Tab_nv_copy[Tab_nv_copy > 0 & Tab_nv_copy < cutoff_variant] <- NA

rows_keep <- NULL
for (i in 1:nrow(Tab_nr_copy)) {
  cat("Working on ", i, "\n")

  x <- Tab_nr_copy[i, ]

  y <- Tab_nv_copy[i, ]

  z <- y

  z[which(x == FALSE)] <- NA

  prop_na <- length(which(is.na(z))) / length(z)

  prop_var <- length(which(z > 0)) / length(z)


  if (prop_na < cutofof_prev_na & prop_var < cutoff_prev_var) {
    rows_keep <- c(rows_keep, i)
  }
}




Tab_nr <- Tab_nr[rows_keep, ]
Tab_nv <- Tab_nv[rows_keep, ]

Tab_nv_copy <- Tab_nv

Tab_nv_copy[Tab_nv_copy < cutoff_variant] <- 0

final_indices <- which(rowSums(Tab_nv_copy, na.rm = TRUE) > 0) %>% as.vector()

Tab_nr <- Tab_nr[final_indices, ]
Tab_nv <- Tab_nv[final_indices, ]

moutfile_nr <- "Tab_NR_somatic_snps_filtered.tsv"
moutfile_nv <- "Tab_NV_somatic_snps_filtered.tsv"

if (mchrs[1] == "chr1") {
  # Write matrices
  mheader <- c("", colnames(Tab_nr))

  write.table(x = mheader %>% t(), file = moutfile_nr, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


  mheader <- c("", colnames(Tab_nv))

  write.table(x = mheader %>% t(), file = moutfile_nv, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
}

write.table(x = Tab_nr, file = moutfile_nr, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)

write.table(x = Tab_nv, file = moutfile_nv, append = TRUE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = FALSE)
