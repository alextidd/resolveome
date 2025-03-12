library(dplyr)
library(magrittr)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

pileup_file <- args[1]

df_pileup <- data.table::fread(pileup_file, sep = "\t", quote = "", header = FALSE) %>%
  as.data.frame()

if (nrow(df_pileup) > 0) {
  df_pileup <- df_pileup %>% dplyr::rename(.data = ., CHROM = V1, POS = V2, REF = V3, NREADS = V4, READ_BASES = V5, BQ = V6, MQ = V7, BPSTART = V8)


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
        ntalt %>% gsub(pattern = "\\|.*", replacement = "")

        indices_top_alt_reverse_a <- indices_norm[(x[indices_norm]) %>% grep(pattern = tolower(ntalt %>% gsub(pattern = "\\|.*", replacement = "")))]

        indices_top_alt_forward_a <- indices_norm[(x[indices_norm]) %>% grep(pattern = toupper(ntalt %>% gsub(pattern = "\\|.*", replacement = "")))]

        indices_top_alt_reverse_b <- indices_norm[(x[indices_norm]) %>% grep(pattern = tolower(ntalt %>% gsub(pattern = ".*\\|", replacement = "")))]

        indices_top_alt_forward_b <- indices_norm[(x[indices_norm]) %>% grep(pattern = toupper(ntalt %>% gsub(pattern = ".*\\|", replacement = "")))]


        # Here mix information
        indices_top_alt_reverse <- c(indices_top_alt_reverse_a, "|", indices_top_alt_reverse_b)

        indices_top_alt_forward <- c(indices_top_alt_forward_a, "|", indices_top_alt_forward_b)
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

  vec_vaf <- NULL
  vec_nvariant <- NULL
  vec_supref <- NULL
  vec_depth <- NULL
  vec_chrom <- NULL
  vec_pos <- NULL
  vec_ref <- NULL
  vec_alt <- NULL

  for (i in 1:nrow(df_pileup)) {
    cat("Working on ", i, " out of ", nrow(df_pileup), "\n")

    malt <- res_decomposition[[i]][[2]]

    x <- df_pileup[i, ]

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

    tograb_ref <- res_decomposition[[i]][[5]]


    if (length(grep(malt, pattern = "\\|")) > 0) {
      vec_chrom <- c(vec_chrom, rep(df_pileup$CHROM[i], 2))
      vec_pos <- c(vec_pos, rep(df_pileup$POS[i], 2))
      vec_ref <- c(vec_ref, rep(df_pileup$REF[i], 2))
      toadd <- c(malt %>% gsub(pattern = "\\|.*", replacement = ""), malt %>% gsub(pattern = ".*\\|", replacement = ""))
      vec_alt <- c(vec_alt, toadd)


      # Now extract information and place it twice
      mforward <- res_decomposition[[i]][[3]]
      mreverse <- res_decomposition[[i]][[4]]


      # Locate the position of the | symbol
      loc_symbol_forward <- which(mforward == "|")

      loc_symbol_reverse <- which(mreverse == "|")

      # Now here check the position
      if (loc_symbol_forward == 1) {
        # Here is the case where there are not positions for first allele
        indices_forward_first_allele <- NA

        indices_forward_second_allele <- mforward %>%
          grep(pattern = "\\|", invert = TRUE, value = TRUE) %>%
          as.numeric()
      } else if (loc_symbol_forward == length(mforward)) {
        # Here is the case where there are not positions for second allele
        indices_forward_first_allele <- mforward %>%
          grep(pattern = "\\|", invert = TRUE, value = TRUE) %>%
          as.numeric()

        indices_forward_second_allele <- NA
      } else {
        # Here we have positions for both alleles
        indices_forward_first_allele <- mforward[1:(loc_symbol_forward - 1)] %>% as.numeric()

        indices_forward_second_allele <- mforward[(loc_symbol_forward + 1):length(mforward)] %>% as.numeric()
      }

      ### Reverse

      if (loc_symbol_reverse == 1) {
        # Here is the case where there are not positions for first allele
        indices_reverse_first_allele <- NA

        indices_reverse_second_allele <- mreverse %>%
          grep(pattern = "\\|", invert = TRUE, value = TRUE) %>%
          as.numeric()
      } else if (loc_symbol_reverse == length(mreverse)) {
        # Here is the case where there are not positions for second allele
        indices_reverse_first_allele <- mreverse %>%
          grep(pattern = "\\|", invert = TRUE, value = TRUE) %>%
          as.numeric()

        indices_reverse_second_allele <- NA
      } else {
        # Here we have positions for both alleles
        indices_reverse_first_allele <- mreverse[1:(loc_symbol_reverse - 1)] %>% as.numeric()

        indices_reverse_second_allele <- mreverse[(loc_symbol_reverse + 1):length(mreverse)] %>% as.numeric()
      }

      # Grab  variants for first allele

      tograb_variant <- c(indices_forward_first_allele, indices_reverse_first_allele)


      indices_variant_mq30_bq30 <- intersect(tograb_variant[which(all_bq[c(tograb_variant)] >= 30)], tograb_variant[which(all_mq[c(tograb_variant)] >= 30)])

      indices_ref_mq30_bq30 <- intersect(tograb_ref[which(all_bq[c(tograb_ref)] >= 30)], tograb_ref[which(all_mq[c(tograb_ref)] >= 30)])

      # Append to type
      vec_depth <- c(vec_depth, length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30))
      vec_nvariant <- c(vec_nvariant, length(indices_variant_mq30_bq30))
      vec_vaf <- c(vec_vaf, length(indices_variant_mq30_bq30) / (length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30)))
      vec_supref <- c(vec_supref, length(indices_ref_mq30_bq30))


      # Grab vriants for second allele

      tograb_variant <- c(indices_forward_second_allele, indices_reverse_second_allele)


      indices_variant_mq30_bq30 <- intersect(tograb_variant[which(all_bq[c(tograb_variant)] >= 30)], tograb_variant[which(all_mq[c(tograb_variant)] >= 30)])

      indices_ref_mq30_bq30 <- intersect(tograb_ref[which(all_bq[c(tograb_ref)] >= 30)], tograb_ref[which(all_mq[c(tograb_ref)] >= 30)])

      # Append to type
      vec_depth <- c(vec_depth, length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30))
      vec_nvariant <- c(vec_nvariant, length(indices_variant_mq30_bq30))
      vec_vaf <- c(vec_vaf, length(indices_variant_mq30_bq30) / (length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30)))
      vec_supref <- c(vec_supref, length(indices_ref_mq30_bq30))
    } else {
      tograb_variant <- c(res_decomposition[[i]][[3]], res_decomposition[[i]][[4]])


      indices_variant_mq30_bq30 <- intersect(tograb_variant[which(all_bq[c(tograb_variant)] >= 30)], tograb_variant[which(all_mq[c(tograb_variant)] >= 30)])

      indices_ref_mq30_bq30 <- intersect(tograb_ref[which(all_bq[c(tograb_ref)] >= 30)], tograb_ref[which(all_mq[c(tograb_ref)] >= 30)])


      vec_depth <- c(vec_depth, length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30))
      vec_nvariant <- c(vec_nvariant, length(indices_variant_mq30_bq30))
      vec_supref <- c(vec_supref, length(indices_ref_mq30_bq30))

      vec_vaf <- c(vec_vaf, length(indices_variant_mq30_bq30) / (length(indices_variant_mq30_bq30) + length(indices_ref_mq30_bq30)))

      vec_chrom <- c(vec_chrom, df_pileup$CHROM[i])
      vec_pos <- c(vec_pos, df_pileup$POS[i])
      vec_ref <- c(vec_ref, df_pileup$REF[i])
      vec_alt <- c(vec_alt, malt)
    }
  }


  ## Create final structure
  df_temp <- data.frame(
    CHROM = vec_chrom,
    POS = vec_pos,
    REF = vec_ref,
    ALT = vec_alt,
    DEPTH = vec_depth,
    NSUPVAR = vec_nvariant,
    NSUPREF = vec_supref,
    VAF = vec_vaf
  )

  write.table(df_temp, file = "res.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
} else {
  df_temp <- data.frame(
    CHROM = NA,
    POS = NA,
    REF = NA,
    ALT = NA,
    DEPTH = NA,
    NSUPVAR = NA,
    NSUPREF = NA,
    VAF = NA
  )

  write.table(df_temp, file = "res.tsv", append = FALSE, quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
}
