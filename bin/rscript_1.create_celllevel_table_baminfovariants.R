library(dplyr)
library(magrittr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (args[1] == "1") {
  pileup_file <- args[2]

  df_pileup <- data.table::fread(pileup_file, sep = "\t", quote = "") %>%
    as.data.frame()

  if (nrow(df_pileup) > 0) {
    df_pileup <- df_pileup %>%
      dplyr::rename(.data = ., CHROM = V1, POS = V2, REF = V3, NREADS = V4, READ_BASES = V5, BQ = V6, MQ = V7, BPSTART = V8, READS = V9, FLAGS = V10, AS = V11)


    # Create function to create segments
    range_vectors <- function(x, y) {
      return(seq(x, y, 1))
    }


    decompose_read_bases <- function(string = NULL) {
      x <- string %>%
        gsub(pattern = "\\\\", replacement = "|") %>%
        strsplit(split = "") %>%
        unlist()

      pos_dol <- x %>% grep(pattern = "\\$")

      # Locate hat positions
      hat_loc <- x %>% grep(pattern = "\\^")
      res_hat <- mapply(range_vectors, hat_loc, hat_loc + 2)
      num_hat <- ncol(res_hat)
      pos_hat <- res_hat %>%
        unlist() %>%
        as.vector()

      # Pos pure dol are end of read matches and the symbol we can omit
      pos_pure_dol <- pos_dol[which(!(pos_dol %in% pos_hat))]

      minus_plus_loc <- x %>% grep(pattern = "[-+]")
      minus_plus_loc <- minus_plus_loc[which(!(minus_plus_loc %in% pos_hat))]

      # Evaluate if here are n


      if (length(minus_plus_loc) > 0) {
        # Here evaluate which positions +2 are numeric and which ones only have one symbol
        minus_plus_loc_two <- minus_plus_loc[which(sapply(x[minus_plus_loc + 2], function(x) x %in% 0:9)) %>% as.vector()]

        minus_plus_loc_one <- minus_plus_loc[which(!(sapply(x[minus_plus_loc + 2], function(x) x %in% 0:9))) %>% as.vector()]


        minus_plus_num_bases_add_two <- paste0(x[minus_plus_loc_two + 1], x[minus_plus_loc_two + 2]) %>% as.numeric()

        minus_plus_num_bases_add_one <- x[minus_plus_loc_one + 1] %>% as.numeric()


        minus_plus_loc <- c(minus_plus_loc_one, minus_plus_loc_two)

        minus_plus_num_bases_add <- c(minus_plus_num_bases_add_one + 1, minus_plus_num_bases_add_two + 2)

        res_mp <- mapply(range_vectors, minus_plus_loc - 1, minus_plus_loc + minus_plus_num_bases_add)


        if (class(res_mp)[1] == "list") {
          num_mp <- length(res_mp)
          pos_mp <- res_mp %>%
            unlist() %>%
            as.vector()
        } else {
          num_mp <- ncol(res_mp)
          pos_mp <- res_mp %>% as.vector()
        }
      } else {
        num_mp <- 0
        pos_mp <- NULL
      }


      star_loc <- x %>% grep(pattern = "\\*")
      pos_star <- star_loc[which(!(star_loc %in% pos_hat))] %>%
        unlist() %>%
        as.vector()
      num_star <- length(pos_star)

      pos_pec <- c(pos_mp, pos_hat, pos_star, pos_pure_dol)

      num_norm <- length(x) - length(pos_pec)

      num_tot <- c(num_norm, num_mp, num_hat, num_star) %>% sum()

      # Now get the indices for the desired allele
      indices_norm <- which(!(seq(1, length(x), 1) %in% pos_pec))

      # Append the indices of the hat
      indices_norm <- c(indices_norm, hat_loc + 2)

      sorted_alt <- x[indices_norm] %>%
        grep(pattern = "[\\.,]", invert = TRUE, value = TRUE) %>%
        toupper() %>%
        table() %>%
        sort(decreasing = TRUE)

      if (length(sorted_alt) == 0) {
        ntalt <- "REF"
      } else if (length(sorted_alt) > 1) {
        if (sorted_alt[1] > sorted_alt[2]) {
          ntalt <- names(sorted_alt[1])
        } else {
          ntalt <- paste0(names(sorted_alt[1]), "|", names(sorted_alt[2]))
        }
      } else {
        ntalt <- names(sorted_alt[1])
      }


      if (ntalt != "REF") {
        if (length(ntalt %>% grep(pattern = "\\|")) > 0) {
          indices_top_alt_reverse <- "MorethanOneAlt"

          indices_top_alt_forward <- "MorethanOneAlt"
        } else {
          # Here distinguish between upper and lower case so we can get strand information
          indices_top_alt_reverse <- indices_norm[(x[indices_norm]) %>% grep(pattern = tolower(ntalt))]

          indices_top_alt_forward <- indices_norm[(x[indices_norm]) %>% grep(pattern = toupper(ntalt))]
        }
      } else {
        indices_top_alt_reverse <- "REF"

        indices_top_alt_forward <- "REF"
      }


      # Now get the reference alleles
      indices_ref <- indices_norm[toupper(x[indices_norm]) %>% grep(pattern = "[\\.,]")]

      return(
        list(
          num_infer_reads = num_tot,
          name_top_alt_al = ntalt,
          indices_top_alt_forward = indices_top_alt_forward,
          indices_top_alt_reverse = indices_top_alt_reverse,
          indices_ref = indices_ref
        )
      )
    }


    res_decomposition <- lapply(df_pileup$READ_BASES, decompose_read_bases)

    num_reads <- res_decomposition %>%
      lapply(function(x) x[1]) %>%
      unlist() %>%
      as.vector()

    table(num_reads == df_pileup$NUM_READS_ALL)

    ## Now process the read indices for the most common alternate allele
    ## We will focus in the specific metrics they detailed in the articles

    vec_median_focal_as <- NULL

    vec_num_fragments_mq40_bq30 <- NULL

    vec_num_forward_mq40_bq30 <- NULL
    vec_sd_forward_mq40_bq30 <- NULL
    vec_mad_forward_mq40_bq30 <- NULL
    vec_prop_forward_mq40_bq30 <- NULL
    vec_num_reverse_mq40_bq30 <- NULL
    vec_sd_reverse_mq40_bq30 <- NULL
    vec_mad_reverse_mq40_bq30 <- NULL
    vec_prop_reverse_mq40_bq30 <- NULL


    vec_depth_mq30_bq30 <- NULL
    vec_vaf_mq30_qb30 <- NULL
    vec_num_var_supporting_mq30_bq30 <- NULL
    vec_num_ref_supporting_mq30_bq30 <- NULL

    vec_gt <- NULL


    vec_names_reads_flags <- NULL
    vec_value_reads_flags <- NULL

    threshold_bp <- (150 * 0.15)

    vec_num_raw_var <- NULL
    vec_num_raw_ref <- NULL

    for (i in 1:nrow(df_pileup)) {
      cat("Working on ", i, " out of ", nrow(df_pileup), "\n")

      x <- df_pileup[i, ]

      tograb_forward <- res_decomposition[[i]][[3]]

      tograb_reverse <- res_decomposition[[i]][[4]]

      ### Alignment Score for reads supporting the variant

      focal_as <- (x$AS %>% strsplit(split = ",") %>% unlist() %>% as.numeric())[c(tograb_forward, tograb_reverse)]

      median_focal_as <- median(focal_as)

      ### Number of clipped reads that support the variant
      focal_reads <- (x$READS %>% strsplit(split = ",") %>% unlist())[c(tograb_forward, tograb_reverse)]

      focal_flags <- (x$FLAGS %>% strsplit(split = ",") %>% unlist())[c(tograb_forward, tograb_reverse)]

      vec_value_reads_flags <- c(vec_value_reads_flags, paste0(focal_reads, "_", focal_flags))

      vec_names_reads_flags <- c(vec_names_reads_flags, rep(paste0("Row", i), length(focal_reads)))

      ### Fragment based statistics are calculated using MQ >=40 and BQ >=40 solely
      all_bq <- x$BQ %>%
        strsplit(split = "") %>%
        unlist() %>%
        sapply(FUN = function(x) utf8ToInt(x) - 33) %>%
        as.vector()

      all_mq <- x$MQ %>%
        strsplit(split = "") %>%
        unlist() %>%
        sapply(FUN = function(x) utf8ToInt(x) - 33) %>%
        as.vector()

      ### Get the high quality fragments supporting the variant
      indices_forward_mq40_bq30 <- intersect(tograb_forward[which(all_bq[c(tograb_forward)] >= 30)], tograb_forward[which(all_mq[c(tograb_forward)] >= 40)])

      indices_reverse_mq40_bq30 <- intersect(tograb_reverse[which(all_bq[c(tograb_reverse)] >= 30)], tograb_reverse[which(all_mq[c(tograb_reverse)] >= 40)])

      # Number of high quality fragments supporting the variant
      num_fragments_mq40_bq30 <- length(indices_forward_mq40_bq30) + length(indices_reverse_mq40_bq30)


      ### Calculate the base SD and MAD using all MQ30 BQ30 reads supporting a variant
      all_bp <- (x$BP %>% strsplit(split = ",") %>% unlist() %>% as.numeric())

      all_bp_mq40_bq30_forward <- all_bp[indices_forward_mq40_bq30]

      all_bp_mq40_bq30_reverse <- all_bp[indices_reverse_mq40_bq30]

      size_forward_mq40_bq30 <- length(all_bp_mq40_bq30_forward)

      size_reverse_mq40_bq30 <- length(all_bp_mq40_bq30_reverse)

      prop_under_forward <- length(which(all_bp_mq40_bq30_forward < threshold_bp)) / size_forward_mq40_bq30

      prop_under_reverse <- length(which(all_bp_mq40_bq30_reverse < threshold_bp)) / size_reverse_mq40_bq30


      sd_forward_mq40_bq30 <- sd(all_bp_mq40_bq30_forward)

      mad_forward_mq40_bq30 <- mad(all_bp_mq40_bq30_forward)


      sd_reverse_mq30_bq30 <- sd(all_bp_mq40_bq30_reverse)

      mad_reverse_mq40_bq30 <- mad(all_bp_mq40_bq30_reverse)


      ### Genotype variants using MQ>=30 and BQ>=30 solely
      indices_forward_mq30_bq30 <- intersect(tograb_forward[which(all_bq[c(tograb_forward)] >= 30)], tograb_forward[which(all_mq[c(tograb_forward)] >= 30)])

      indices_reverse_mq30_bq30 <- intersect(tograb_reverse[which(all_bq[c(tograb_reverse)] >= 30)], tograb_reverse[which(all_mq[c(tograb_reverse)] >= 30)])

      indices_variant_mq30_bq30 <- c(indices_forward_mq30_bq30, indices_reverse_mq30_bq30)


      tograb_ref <- res_decomposition[[i]][[5]]


      indices_ref_mq30_bq30 <- intersect(tograb_ref[which(all_bq[c(tograb_ref)] >= 30)], tograb_ref[which(all_mq[c(tograb_ref)] >= 30)])


      depth_mq30_bq30 <- length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30)

      vaf_mq30_bq30 <- length(indices_variant_mq30_bq30) / depth_mq30_bq30

      num_var_supporting_mq30_bq30 <- length(indices_variant_mq30_bq30)

      ### Fill vectors with results
      vec_median_focal_as <- c(vec_median_focal_as, median_focal_as)

      vec_num_fragments_mq40_bq30 <- c(vec_num_fragments_mq40_bq30, num_fragments_mq40_bq30)

      vec_num_forward_mq40_bq30 <- c(vec_num_forward_mq40_bq30, size_forward_mq40_bq30)
      vec_sd_forward_mq40_bq30 <- c(vec_sd_forward_mq40_bq30, sd_forward_mq40_bq30)
      vec_mad_forward_mq40_bq30 <- c(vec_mad_forward_mq40_bq30, mad_forward_mq40_bq30)
      vec_prop_forward_mq40_bq30 <- c(vec_prop_forward_mq40_bq30, prop_under_forward)

      vec_num_reverse_mq40_bq30 <- c(vec_num_reverse_mq40_bq30, size_reverse_mq40_bq30)
      vec_sd_reverse_mq40_bq30 <- c(vec_sd_reverse_mq40_bq30, sd_reverse_mq30_bq30)
      vec_mad_reverse_mq40_bq30 <- c(vec_mad_reverse_mq40_bq30, mad_reverse_mq40_bq30)
      vec_prop_reverse_mq40_bq30 <- c(vec_prop_reverse_mq40_bq30, prop_under_reverse)

      vec_depth_mq30_bq30 <- c(vec_depth_mq30_bq30, depth_mq30_bq30)
      vec_vaf_mq30_qb30 <- c(vec_vaf_mq30_qb30, vaf_mq30_bq30)
      vec_num_var_supporting_mq30_bq30 <- c(vec_num_var_supporting_mq30_bq30, num_var_supporting_mq30_bq30)
      vec_num_ref_supporting_mq30_bq30 <- c(vec_num_ref_supporting_mq30_bq30, length(indices_ref_mq30_bq30))
      vec_gt <- c(vec_gt, res_decomposition[[i]][[2]])

      vec_num_raw_var <- c(vec_num_raw_var, length(tograb_forward) + length(tograb_reverse))
      vec_num_raw_ref <- c(vec_num_raw_ref, length(tograb_ref))
    }


    df_meta <- df_pileup[, 1:3]
    rm(df_pileup)
    rm(res_decomposition)
    gc()

    df_read_flags <- data.frame(
      ID_FLAG = vec_value_reads_flags,
      IdRow = vec_names_reads_flags
    )


    rm(vec_value_reads_flags)
    rm(vec_names_reads_flags)


    ### Create data frame with results and print table
    res <- df_meta %>%
      cbind(
        data.frame(
          MEDIAN_AS_VARIANT_READS = vec_median_focal_as,
          NUM_FRAGMENTS_MQ40_BQ30_FR = vec_num_fragments_mq40_bq30,
          NUM_FRAGMENTS_MQ40_BQ30_F = vec_num_forward_mq40_bq30,
          NUM_FRAGMENTS_MQ40_BQ30_R = vec_num_reverse_mq40_bq30,
          PROP_FRAGMENTS_BPSTART_UNDER_MQ40_BQ30_F = vec_prop_forward_mq40_bq30,
          PROP_FRAGMENTS_BPSTART_UNDER_MQ40_BQ30_R = vec_prop_reverse_mq40_bq30,
          SD_BPSTART_FRAGMENTS_MQ40_BQ30_F = vec_sd_forward_mq40_bq30,
          SD_BPSTART_FRAGMENTS_MQ40_BQ30_R = vec_sd_reverse_mq40_bq30,
          MAD_BPSTART_RAGMENTS_MQ40_BQ30_F = vec_mad_forward_mq40_bq30,
          MAD_BPSTART_RAGMENTS_MQ40_BQ30_R = vec_mad_reverse_mq40_bq30,
          DEPTH_MQ30_BQ30 = vec_depth_mq30_bq30,
          VAF_MQ30_BQ30 = vec_vaf_mq30_qb30,
          NUM_READS_REF_READS_MQ30_BQ30 = vec_num_ref_supporting_mq30_bq30,
          NUM_READS_VARIANT_READS_MQ30_BQ30 = vec_num_var_supporting_mq30_bq30,
          NUM_READS_REF_READS_MQ0_BQ0 = vec_num_raw_ref,
          NUM_READS_VARIANT_READS_MQ0_BQ0 = vec_num_raw_var,
          ALT_GT_INFERRED = vec_gt
        )
      ) %>%
      dplyr::mutate(.data = ., IdRow = paste0("Row", 1:nrow(df_meta)))

    rm(vec_median_focal_as)
    rm(vec_num_fragments_mq40_bq30)
    rm(vec_num_forward_mq40_bq30)
    rm(vec_num_reverse_mq40_bq30)
    rm(vec_prop_forward_mq40_bq30)
    rm(vec_prop_reverse_mq40_bq30)
    rm(vec_sd_forward_mq40_bq30)
    rm(vec_sd_reverse_mq40_bq30)
    rm(vec_mad_forward_mq40_bq30)
    rm(vec_mad_reverse_mq40_bq30)
    rm(vec_depth_mq30_bq30)
    rm(vec_vaf_mq30_qb30)
    rm(vec_num_ref_supporting_mq30_bq30)
    rm(vec_num_var_supporting_mq30_bq30)
    rm(vec_gt)
    rm(vec_num_raw_var)
    rm(vec_num_raw_ref)

    gc()


    ### Write table
    write.table(res, file = "res.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

    write.table(df_read_flags, file = "res_flags.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  }
} else {
  df_name_flag_cigar <- data.table::fread("df_found", sep = "\t", header = FALSE) %>%
    as.data.frame() %>%
    dplyr::rename(.data = ., ID = V1, FLAG = V2, CIGAR = V3) %>%
    dplyr::mutate(.data = ., ID_FLAG = paste0(ID, "_", FLAG))

  df_name_flag_cigar$BooleanClipped <- 0
  df_name_flag_cigar$BooleanClipped[df_name_flag_cigar$CIGAR %>% grep(pattern = "H|S")] <- 1

  df_read_flags <- data.table::fread("res_flags.tsv", sep = "\t", header = FALSE) %>%
    as.data.frame() %>%
    dplyr::rename(.data = ., ID_FLAG = V1, IdRow = V2)

  df_read_flags$BooleanClipped <- match(df_read_flags$ID_FLAG, df_name_flag_cigar$ID_FLAG) %>% df_name_flag_cigar$BooleanClipped[.]

  rm(df_name_flag_cigar)
  gc()

  df_read_flags$Dummy <- 1

  df_sum <- aggregate(BooleanClipped ~ IdRow, df_read_flags, sum)

  df_tot <- aggregate(Dummy ~ IdRow, df_read_flags, sum)

  merged <- merge(df_tot, df_sum, by = "IdRow", all = TRUE) %>%
    dplyr::mutate(.data = ., PropClippedSupportVariant = BooleanClipped / Dummy) %>%
    dplyr::select(.data = ., -c("BooleanClipped", "Dummy"))

  rm(df_sum)
  rm(df_tot)
  gc()

  res <- read.table(file = "res.tsv", header = TRUE, sep = "\t", quote = "", comment.char = "")


  # Merge structure into final data frame
  res <- merge(merged, res, by = "IdRow")

  write.table(res, file = "res_end.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
